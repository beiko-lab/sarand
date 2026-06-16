"""Parsing of assembly-graph node paths and Bandage/BLAST alignment output."""
from __future__ import annotations

import collections
import csv
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

from sarand.external.bandage import BandageResult
from sarand.util.file import try_dump_to_disk
from sarand.util.logger import LOG


def reverse_sign(sign: str) -> str:
    """Reverse a strand sign (+/-)."""
    if sign == "-":
        return "+"
    elif sign == "+":
        return "-"
    else:
        LOG.error("ERROR: ivalid sign!")
        sys.exit(1)


def find_node_name(node: str) -> str:
    """
    To remove specific characters and return the rest except the last character
    as the node name
    """
    return re.sub("[]{}[]", "", node)[:-1]


def find_node_name_orient(node: str) -> str:
    """
    To remove specific characters and return the rest as the name+orient of the node
    """
    return re.sub("[]{}[]", "", node)


def exist_in_path(path: List[str], mynode: str) -> int:
    """
    To check if a given node exists in the path
    Parameters:
        path:	a list of nodes
        mynode:	the node to check if it is in the path
    Return:
        the index of mynode in the path if present; -1 otherwise
    """
    for i, node in enumerate(path):
        if find_node_name_orient(node) == mynode:
            return i
    return -1


def extract_nodes_in_path(path: str) -> Tuple[List[str], List[str], int, int]:
    """
    Parameters:
        path:	a list of nodes with -/+ tail and comma separated : e.g., '(1363) 69-, 2193+ (1786)'
    Return:
        node_list:	list of node numbers -> e.g., [69, 2193]
        orientation_list: list of orientation of nodes -> e.g., [-, +]
        start_pos:	where in the first node the sequence has started -> e.g., 1363
        end_pos:	where in the last node the sequence ended  -> e.g., 1786
    """
    start_pos = 0
    end_pos = 0
    if path.startswith("("):
        index = path.find(")")
        start_pos = int(path[1:index])
    if path.endswith(")"):
        index = path.rfind("(")
        end_pos = int(path[index + 1: -1])
    # Remove text between ()
    path = (re.sub(r"\((.*?)\)", "", path)).strip()
    node_list = []
    orientation_list = []
    nodes = path.split(",")
    for node in nodes:
        if "-" in node:
            orientation_list.append("-")
        else:
            orientation_list.append("+")
        node = re.sub("[+-]", "", node.split()[0])
        node_list.append(node)
    return node_list, orientation_list, start_pos, end_pos


def read_path_info_from_align_file(
        align_file: str | Path, threshold: float = 95
) -> Tuple[bool, List[Dict[str, Any]]]:
    """Parse a single-AMR Bandage alignment TSV into per-path node/position info.

    Parameters:
        align_file: the Bandage querypaths ``.tsv`` output.
        threshold: minimum coverage and identity percentage to keep a path.
    Return:
        (found, paths_info) where ``found`` is True if any path passed the
        threshold and ``paths_info`` is a list of node/orientation/position dicts.
    """
    paths_info = []
    found = False
    with open(align_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        # skip the header
        next(reader)
        for row in reader:
            coverage = float(re.sub("[%]", "", row[3]))
            identity = float(re.sub("[%]", "", row[5]))
            if int(coverage) >= threshold and int(identity) >= threshold:
                found = True
                cell_info = row[1].strip()
                nodes, orientation_list, start_pos, end_pos = extract_nodes_in_path(
                    cell_info
                )
                path_info = {
                    "nodes": nodes,
                    "orientations": orientation_list,
                    "start_pos": start_pos,
                    "end_pos": end_pos,
                }
                paths_info.append(path_info)
    if not found:
        LOG.error("No path info was found in " + align_file)
    return found, paths_info


def read_path_info_from_align_file_with_multiple_amrs(
        output_name: Path,
        ga: List[BandageResult],
        threshold: float = 99,
        debug: bool = False
) -> Dict[str, List[Dict[str, Any]]]:
    """Group Bandage results by AMR into per-path node/position info.

    Parameters:
        output_name: directory used for optional debug output.
        ga: the parsed Bandage results to filter and group.
        threshold: minimum coverage and identity percentage to keep a path.
        debug: if True, write the coverage/identity of every result to disk.
    Return:
        a mapping of AMR name -> list of node/orientation/position dicts.
    """
    debug_to_write = list()

    paths_info_list = collections.defaultdict(list)
    for result in ga:
        amr_name = result.amr_name
        coverage = result.coverage_pct
        identity = result.identity_pct

        # Append debugging information
        if debug:
            debug_to_write.append({
                'id': result.identity,
                'path': f'({result.path_start}) {result.path} ({result.path_end})',
                'amr': amr_name,
                'coverage': coverage,
                'identity': identity
            })

        if coverage >= threshold and identity >= threshold:
            nodes, orientation_list = result.path_to_sarand
            start_pos = result.path_start
            end_pos = result.path_end
            path_info = {
                "nodes": nodes,
                "orientations": orientation_list,
                "start_pos": start_pos,
                "end_pos": end_pos,
            }
            paths_info_list[amr_name].append(path_info)

    if debug:
        try_dump_to_disk(
            debug_to_write,
            output_name / 'aligner_tool_coverage_identity.json'
        )
    return paths_info_list
