"""Parsing of assembly-graph node paths and Bandage/BLAST alignment output."""
from __future__ import annotations

import collections
from pathlib import Path
from typing import Any, Dict, List

from sarand.external.bandage import BandageResult
from sarand.util.file import try_dump_to_disk


def read_path_info_from_align_file_with_multiple_targets(
        output_name: Path,
        ga: List[BandageResult],
        min_target_identity: float = 95,
        min_target_coverage: float = 95,
        debug: bool = False
) -> Dict[str, List[Dict[str, Any]]]:
    """Group Bandage results by target into per-path node/position info.

    Parameters:
        output_name: directory used for optional debug output.
        ga: the parsed Bandage results to filter and group.
        threshold: minimum coverage and identity percentage to keep a path.
        debug: if True, write the coverage/identity of every result to disk.
    Return:
        a mapping of target name -> list of node/orientation/position dicts.
    """
    debug_to_write = list()

    paths_info_list = collections.defaultdict(list)
    for result in ga:
        target_name = result.target_name
        coverage = result.coverage_pct
        identity = result.identity_pct

        # Append debugging information
        if debug:
            debug_to_write.append({
                'id': result.identity,
                'path': f'({result.path_start}) {result.path} ({result.path_end})',
                'target': target_name,
                'coverage': coverage,
                'identity': identity
            })

        if coverage >= min_target_coverage and identity >= min_target_identity:
            nodes, orientation_list = result.path_to_sarand
            start_pos = result.path_start
            end_pos = result.path_end
            path_info = {
                "nodes": nodes,
                "orientations": orientation_list,
                "start_pos": start_pos,
                "end_pos": end_pos,
            }
            paths_info_list[target_name].append(path_info)

    if debug:
        try_dump_to_disk(
            debug_to_write,
            output_name / 'aligner_tool_coverage_identity.json'
        )
    return paths_info_list
