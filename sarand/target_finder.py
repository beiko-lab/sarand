"""Stage 1 of the graph pipeline: locate target genes in the assembly graph.

Target genes are aligned to the graph with Bandage+BLAST (grouped and run in a
multiprocessing pool), hits with overlapping paths are de-duplicated, and the
unique hits are written to disk.
"""
from __future__ import annotations

import collections
import datetime
from functools import partial
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from sarand.config import (
    TARGET_ALIGN_DIR,
    TARGET_DIR_NAME,
    TARGET_OVERLAP_FILE,
    TARGET_SEQ_DIR,
)
from sarand.external.bandage import Bandage
from sarand.util.file import try_dump_to_disk
from sarand.util.graph_path import read_path_info_from_align_file_with_multiple_targets
from sarand.util.logger import LOG
from sarand.util.naming import (
    target_name_from_comment,
    extract_name_from_file_name,
    restricted_target_name_from_modified_name,
)
from sarand.util.sequence import create_fasta_file


def target_path_overlap(
        found_target_paths: List[List[Dict[str, Any]]],
        new_paths: List[Dict[str, Any]],
        new_target_len: int,
        overlap_percent: int = 95,
) -> Tuple[bool, Optional[List[int]]]:
    """
    To check if all paths found for the new target seq overlap significantly (greater/equal
     than/to overlap percent) with the already found paths for other targets 
    Parameters:
         found_target_paths:  	the paths already found for target genes
        new_paths: 			the paths found for the new target gene
        overlap_percent:	the threshold for overlap
    Return:
        False only if every paths in new_paths have overlap with at least one path in found_target_paths
        True if we can find at least one path that is unique and not available in found_target_paths
        Also, it returns the list of indeces from found_target_paths that had overlap with a path in new_paths
    """
    id_list = []
    for new_path in new_paths:
        found = False
        for i, paths in enumerate(found_target_paths):
            for path in paths:
                if (
                        path["nodes"] == new_path["nodes"]
                        and path["orientations"] == new_path["orientations"]
                ):
                    # for now we just check overlaps when they are in the same node(s)
                    diff_length = max(
                        path["start_pos"] - new_path["start_pos"], 0
                    ) + max(new_path["end_pos"] - path["end_pos"], 0)
                    percent = (1 - (float(diff_length) / (new_target_len))) * 100
                    if percent >= overlap_percent:
                        found = True
                        if i not in id_list:
                            id_list.append(i)
                        break
            if found:
                break
    if len(id_list) == len(new_paths):
        return True, id_list
    return False, None


def align_targets_to_graph(
        gfa_file: Path,
        output_dir: Path,
        min_target_identity: float,
        min_target_coverage: float,
        target_object: Tuple[Path, List[str]],
        keep_files: bool,
        debug: bool,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Align a group of target sequences to the assembly graph with Bandage+BLAST and
    return, per target, the graph paths whose coverage and identity meet the
    thresholds.
    Parameters:
        gfa_file: the address of the assembly graph
        output_dir: the address of the output directory
        threshold: the threshold for coverage and identity
        min_target_identity: float,
        min_target_coverage: float,
        target_object: the concatenated AMR file and the list of AMR file names
        keep_files: True if intermediate files should be kept, False otherwise.
        debug: True if additional debug files should be created, False otherwise.
    """
    cat_file, target_files = target_object
    target_names = [extract_name_from_file_name(e) for e in target_files]
    LOG.debug(
        'Checking if target gene "' + str(target_names) + '" can be found in the assembly graph...'
    )
    output_name = Path(output_dir) / (
        extract_name_from_file_name(cat_file)
        + "_align_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    )

    # Run Bandage+BLAST
    aligner_path = output_name if keep_files else None
    aligner = Bandage.run_for_sarand(
        gfa=gfa_file,
        reads=cat_file,
        min_target_identity=min_target_identity,
        min_target_coverage=min_target_coverage,
        out_dir=aligner_path,
    )
    
    paths_info_list = read_path_info_from_align_file_with_multiple_targets(
        output_name=output_dir,
        ga=aligner.results,
        min_target_identity=min_target_identity,
        min_target_coverage=min_target_coverage,
        debug=debug,
    )

    return paths_info_list


def find_target_group_in_graph(
        gfa_file: Path,
        align_dir: Path,
        output_dir: Path,
        min_target_identity: float,
        min_target_coverage: float,
        keep_files: bool,
        core_num: int,
        debug: bool,
        target_object: Tuple[int, List[Tuple[str, str]]],
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Write a group of target sequences into a single file and align it against the graph.
    Parameters:
        gfa_file: the file containing the assembly graph
        align_dir: the directory for storing alignment info
        output_dir: the directory to store the list of AMRs in a single file
        min_target_identity: minimum target identity
        min_target_coverage: minimum target coverage
        keep_files: True if intermediate files should be kept, False otherwise.
        core_num: the number of cores used
        debug: True if additional debug files should be created, False otherwise.
        target_object: the list of AMRs and their ids
    Return:
        the alignment info for AMRs
    """
    g_id, target_group = target_object
    # read info of the group into a single file
    cat_file = Path(output_dir) / TARGET_DIR_NAME / ("target_group_" + str(g_id) + ".fasta")
    file_group = []
    with open(cat_file, "w") as writer:
        for target_info in target_group:
            target_seq, target_title = target_info
            writer.write(target_title)
            writer.write(target_seq)
            target_name1 = target_name_from_comment(target_title)
            target_file_name = restricted_target_name_from_modified_name(target_name1)
            file_group.append(target_file_name + ".fasta")

    # Run Bandage+BLAST
    p_find_target_align = align_targets_to_graph(
        gfa_file=gfa_file,
        output_dir=Path(align_dir),
        min_target_identity=min_target_identity,
        min_target_coverage=min_target_coverage,
        target_object=(cat_file, file_group),
        keep_files=keep_files,
        debug=debug,
    )
    if debug:
        try_dump_to_disk(p_find_target_align, Path(align_dir) / 'debug_p_find_target_align.json')

    # Remove temporary target file
    if cat_file.is_file():
        cat_file.unlink()
    return p_find_target_align


def collect_unique_targets(
        target_seq_title_list: List[Tuple[str, str]],
        target_group_id: Dict[str, int],
        paths_info_group_list: List[Dict[str, List[Dict[str, Any]]]],
) -> Tuple[List[str], List[Dict[str, Any]], List[List[Dict[str, Any]]]]:
    """
    Process the per-group alignment results, keeping targets with unique graph paths
    and group targets whose paths overlap.
    """
    unique_target_seqs = []
    unique_target_infos = []
    unique_target_paths = []
    for target_object in target_seq_title_list:
        target_name = target_name_from_comment(target_object[1])
        id = target_group_id[target_name]
        restricted_target_name = restricted_target_name_from_modified_name(target_name)
        if restricted_target_name in paths_info_group_list[id]:
            LOG.debug(target_name + " was found!")
            path_info = paths_info_group_list[id][restricted_target_name]
            overlap, target_ids = target_path_overlap(
                unique_target_paths, path_info, len(target_object[0]) - 1
            )
            if not overlap:
                unique_target_seqs.append(target_object[0])
                target_info = {"name": target_object[1], "overlap_list": []}
                unique_target_infos.append(target_info)
                unique_target_paths.append(path_info)
            else:
                if len(target_ids) > 1:
                    LOG.error("A target hit has overlap with more than one group")
                # add this target to the right group of targets all sharing overlaps
                for id in target_ids:
                    if target_name not in unique_target_infos[id]["overlap_list"]:
                        unique_target_infos[id]["overlap_list"].append(target_name)

    return unique_target_seqs, unique_target_infos, unique_target_paths


def find_all_targets_in_graph(
        gfa_file: Path,
        output_dir: str,
        target_sequences_file: Path,
        min_target_identity: float,
        min_target_coverage: float,
        core_num: int,
        keep_files: bool,
        debug: bool,
) -> Tuple[List[str], List[List[Dict[str, Any]]]]:
    """
    Align every target gene in ``target_sequences_file`` against the graph and
    return the unique hits and their paths.
    Parameters:
        gfa_file: the address of the assembly graph
        output_dir: the address of the output directory
        target_sequences_file: the address of the file containing the sequence of all AMRs from CARD
        min_target_identity: minimum target identity
        min_target_coverage: minimum target coverage
        core_num: the number of used cores
        keep_files: True if intermediate files should be kept, False otherwise.
        debug: True if additional debug files should be created, False otherwise.
    """
    align_dir = Path(output_dir) / TARGET_DIR_NAME / TARGET_ALIGN_DIR
    align_dir.mkdir(parents=True, exist_ok=True)

    # generate the groups and store the group of each target
    group_num = 5
    target_group_id = collections.defaultdict(list)
    target_file_groups = [[] for i in range(group_num * core_num)]
    target_title = ""
    target_seq_title_list = []
    # Read target sequences one by one
    target_counter = 0
    with open(target_sequences_file) as fp:
        for line in fp:
            if line.startswith(">"):
                target_title = line
                continue
            target_name = target_name_from_comment(target_title[:-1])
            target_seq_title_list.append((line, target_title))
            id = target_counter % (group_num * core_num)
            target_file_groups[id].append((line, target_title))
            target_group_id[target_name] = id
            target_counter += 1

    target_objects = [(i, e) for i, e in enumerate(target_file_groups)]
    # parallel run Bandage+BLAST
    p_find_target = partial(
        find_target_group_in_graph,
        gfa_file,
        Path(align_dir),
        Path(output_dir),
        min_target_identity,
        min_target_coverage,
        keep_files,
        core_num,
        debug,
    )
    with Pool(core_num) as p:
        paths_info_group_list = p.map(p_find_target, target_objects)
        
    unique_target_seqs, unique_target_infos, unique_target_paths = collect_unique_targets(
        target_seq_title_list, target_group_id, paths_info_group_list)

    if debug:
        try_dump_to_disk(
            [
                {'unique_target_seqs': s, 'unique_target_infos': i, 'unique_target_paths': p}
                for s, i, p in zip(unique_target_seqs, unique_target_infos, unique_target_paths)
            ],
            Path(align_dir) / 'debug_get_unique_target_info.json'
        )

    # write the sequence of found target seqs that don't have overlapped paths with
    # others + the list of groups in which all targets have overlapping paths
    unique_target_files = write_found_targets_to_disk(
        output_dir=Path(output_dir),
        unique_target_seqs=unique_target_seqs,
        unique_target_infos=unique_target_infos
    )

    return unique_target_files, unique_target_paths


def write_found_targets_to_disk(
        output_dir: Path,
        unique_target_seqs: List[str],
        unique_target_infos: List[Dict[str, Any]],
) -> List[str]:
    """Write each unique target sequence to its own FASTA and record overlap groups.

    Parameters:
        output_dir: the run output directory.
        unique_target_seqs: the sequences of the de-duplicated target sequence hits.
        unique_target_infos: per-hit metadata (name + list of overlapping target hits).
    Return:
        the list of written per-target FASTA file paths.
    """
    target_dir = output_dir / TARGET_DIR_NAME / TARGET_SEQ_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    overlap_file_name = output_dir / TARGET_DIR_NAME / TARGET_OVERLAP_FILE

    unique_target_files = list()

    with overlap_file_name.open('w') as f:
        for i, seq in enumerate(unique_target_seqs):
            target_name = target_name_from_comment(unique_target_infos[i]["name"])
            restricted_target_name = restricted_target_name_from_modified_name(target_name)
            target_file = create_fasta_file(
                seq,
                str(target_dir.absolute()),
                f'>{unique_target_infos[i]["name"]}',
                restricted_target_name
            )
            unique_target_files.append(target_file)
            f.write(target_name + ":")
            if unique_target_infos[i]["overlap_list"]:
                f.write(", ".join(e for e in unique_target_infos[i]["overlap_list"]))
                f.write("\n")
            else:
                f.write("\n")
    return unique_target_files
