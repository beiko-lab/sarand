"""Stage 1 of the graph pipeline: locate target (AMR) genes in the assembly graph.

Target genes are aligned to the graph with Bandage+BLAST (grouped and run in a
multiprocessing pool), hits with overlapping paths are de-duplicated, and the
unique hits are written to disk.
"""
import collections
import datetime
import os
from functools import partial
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Any, Dict, List

from sarand.config import (
    AMR_ALIGN_DIR,
    AMR_DIR_NAME,
    AMR_OVERLAP_FILE,
    AMR_SEQ_DIR,
)
from sarand.external.bandage import Bandage
from sarand.util.file import try_dump_to_disk
from sarand.util.graph_path import read_path_info_from_align_file_with_multiple_amrs
from sarand.util.logger import LOG
from sarand.util.naming import (
    amr_name_from_comment,
    extract_name_from_file_name,
    restricted_amr_name_from_modified_name,
)
from sarand.util.sequence import create_fasta_file


def amr_path_overlap(found_amr_paths, new_paths, new_amr_len, overlap_percent=95):
    """
    To check if all paths found for the new AMR seq overlap significantly (greater/equal
     than/to overlap percent) with the already found paths for other AMRs
    Parameters:
         found_amr_paths:  	the paths already found for AMR genes
        new_paths: 			the paths found for the new AMR gene
        overlap_percent:	the threshold for overlap
    Return:
        False only if every paths in new_paths have overlap with at least one path in found_amr_paths
        True if we can find at least one path that is unique and not available in found_amr_paths
        Also, it returns the list of indeces from found_amr_paths that had overlap with a path in new_paths
    """
    id_list = []
    for new_path in new_paths:
        found = False
        for i, paths in enumerate(found_amr_paths):
            for path in paths:
                if (
                        path["nodes"] == new_path["nodes"]
                        and path["orientations"] == new_path["orientations"]
                ):
                    # for now we just check overlaps when they are in the same node(s)
                    diff_length = max(
                        path["start_pos"] - new_path["start_pos"], 0
                    ) + max(new_path["end_pos"] - path["end_pos"], 0)
                    percent = (1 - (float(diff_length) / (new_amr_len))) * 100
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


def align_amrs_to_graph(
        gfa_file: Path,
        output_dir: Path,
        threshold: float,
        amr_object,
        keep_files: bool,
        debug: bool,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Align a group of AMR sequences to the assembly graph with Bandage+BLAST and
    return, per AMR, the graph paths whose coverage and identity meet the
    threshold.
    Parameters:
        gfa_file: the address of the assembly graph
        output_dir: the address of the output directory
        threshold: the threshold for coverage and identity
        amr_object: the concatenated AMR file and the list of AMR file names
        keep_files: True if intermediate files should be kept, False otherwise.
        debug: True if additional debug files should be created, False otherwise.
    """
    cat_file, amr_files = amr_object
    amr_names = [extract_name_from_file_name(e) for e in amr_files]
    LOG.debug(
        'Checking if AMR gene "' + str(amr_names) + '" exists in the assembly graph...'
    )
    output_name = os.path.join(
        output_dir,
        extract_name_from_file_name(cat_file)
        + "_align_"
        + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
        "")

    # Run Bandage+BLAST
    aligner_path = output_name if keep_files else None
    aligner = Bandage.run_for_sarand(
        gfa=gfa_file,
        reads=cat_file,
        threshold=threshold,
        out_dir=aligner_path,
    )

    paths_info_list = read_path_info_from_align_file_with_multiple_amrs(
        output_name=output_dir,
        ga=aligner.results,
        threshold=threshold,
        debug=debug,
    )

    return paths_info_list


def find_amr_group_in_graph(
        gfa_file: Path,
        align_dir: Path,
        output_dir: Path,
        amr_threshold: float,
        keep_files: bool,
        core_num: int,
        debug: bool,
        amr_object
):
    """
    Write a group of AMRs into a single file and align it against the graph.
    Parameters:
        gfa_file: the file containing the assembly graph
        align_dir: the directory for storing alignment info
        output_dir: the directory to store the list of AMRs in a single file
        amr_threshold: the threshold for identity and coverage used in alignment
        keep_files: True if intermediate files should be kept, False otherwise.
        core_num: the number of cores used
        debug: True if additional debug files should be created, False otherwise.
        amr_object: the list of AMRs and their ids
    Return:
        the alignment info for AMRs
    """
    g_id, amr_group = amr_object
    # read info of the group into a single file
    cat_file = os.path.join(
        output_dir, AMR_DIR_NAME, "amr_group_" + str(g_id) + ".fasta"
    )
    file_group = []
    with open(cat_file, "w") as writer:
        for amr_info in amr_group:
            amr_seq, amr_title = amr_info
            writer.write(amr_title)
            writer.write(amr_seq)
            amr_name1 = amr_name_from_comment(amr_title)
            amr_file_name = restricted_amr_name_from_modified_name(amr_name1)
            file_group.append(amr_file_name + ".fasta")

    # Run Bandage+BLAST
    p_find_amr_align = align_amrs_to_graph(
        gfa_file=gfa_file,
        output_dir=Path(align_dir),
        threshold=amr_threshold,
        amr_object=(cat_file, file_group),
        keep_files=keep_files,
        debug=debug,
    )
    if debug:
        try_dump_to_disk(p_find_amr_align, Path(align_dir) / 'debug_p_find_amr_align.json')

    # Remove temporary AMR file
    if os.path.isfile(cat_file):
        os.remove(cat_file)
    return p_find_amr_align


def collect_unique_amrs(amr_seq_title_list, amr_group_id, paths_info_group_list):
    """
    Process the per-group alignment results, keeping AMRs with unique graph paths
    and grouping AMRs whose paths overlap.
    """
    unique_amr_seqs = []
    unique_amr_infos = []
    unique_amr_paths = []
    for amr_object in amr_seq_title_list:
        amr_name = amr_name_from_comment(amr_object[1])
        id = amr_group_id[amr_name]
        restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
        if restricted_amr_name in paths_info_group_list[id]:
            LOG.debug(amr_name + " was found!")
            path_info = paths_info_group_list[id][restricted_amr_name]
            overlap, amr_ids = amr_path_overlap(
                unique_amr_paths, path_info, len(amr_object[0]) - 1
            )
            if not overlap:
                unique_amr_seqs.append(amr_object[0])
                amr_info = {"name": amr_object[1], "overlap_list": []}
                unique_amr_infos.append(amr_info)
                unique_amr_paths.append(path_info)
            else:
                if len(amr_ids) > 1:
                    LOG.error("an AMR has overlap with more than one group")
                # add this AMR to the right group of AMRs all having overlaps
                for id in amr_ids:
                    if amr_name not in unique_amr_infos[id]["overlap_list"]:
                        unique_amr_infos[id]["overlap_list"].append(amr_name)

    return unique_amr_seqs, unique_amr_infos, unique_amr_paths


def find_all_amr_in_graph(
        gfa_file: Path,
        output_dir: str,
        amr_sequences_file: Path,
        amr_threshold: float,
        core_num: int,
        keep_files: bool,
        debug: bool,
):
    """
    Align every target gene in ``amr_sequences_file`` against the graph and
    return the unique hits and their paths.
    Parameters:
        gfa_file: the address of the assembly graph
        output_dir: the address of the output directory
        amr_sequences_file: the address of the file containing the sequence of all AMRs from CARD
        amr_threshold: the threshold for coverage and identity
        core_num: the number of used cores
        keep_files: True if intermediate files should be kept, False otherwise.
        debug: True if additional debug files should be created, False otherwise.
    """
    align_dir = os.path.join(output_dir, AMR_DIR_NAME, AMR_ALIGN_DIR)
    os.makedirs(align_dir, exist_ok=True)

    # generate the groups and store the group of each amr
    group_num = 5
    amr_group_id = collections.defaultdict(list)
    amr_file_groups = [[] for i in range(group_num * core_num)]
    amr_title = ""
    amr_seq_title_list = []
    # Read AMR sequences one by one
    amr_counter = 0
    with open(amr_sequences_file) as fp:
        for line in fp:
            if line.startswith(">"):
                amr_title = line
                continue
            amr_name = amr_name_from_comment(amr_title[:-1])
            amr_seq_title_list.append((line, amr_title))
            id = amr_counter % (group_num * core_num)
            amr_file_groups[id].append((line, amr_title))
            amr_group_id[amr_name] = id
            amr_counter += 1

    amr_objects = [(i, e) for i, e in enumerate(amr_file_groups)]
    # parallel run Bandage+BLAST
    p_find_amr = partial(
        find_amr_group_in_graph,
        gfa_file,
        Path(align_dir),
        Path(output_dir),
        amr_threshold,
        keep_files,
        core_num,
        debug,
    )
    with Pool(core_num) as p:
        paths_info_group_list = p.map(p_find_amr, amr_objects)

    unique_amr_seqs, unique_amr_infos, unique_amr_paths = collect_unique_amrs(
        amr_seq_title_list, amr_group_id, paths_info_group_list)

    if debug:
        try_dump_to_disk(
            [
                {'unique_amr_seqs': s, 'unique_amr_infos': i, 'unique_amr_paths': p}
                for s, i, p in zip(unique_amr_seqs, unique_amr_infos, unique_amr_paths)
            ],
            Path(align_dir) / 'debug_get_unique_amr_info.json'
        )

    # write the sequence of found AMRs that don't have overlapped paths with
    # others + the list of groups in which all AMRs have overlapped paths
    unique_amr_files = write_found_amrs_to_disk(
        output_dir=Path(output_dir),
        unique_amr_seqs=unique_amr_seqs,
        unique_amr_infos=unique_amr_infos
    )

    return unique_amr_files, unique_amr_paths


def write_found_amrs_to_disk(output_dir: Path, unique_amr_seqs, unique_amr_infos):
    amr_dir = output_dir / AMR_DIR_NAME / AMR_SEQ_DIR
    os.makedirs(amr_dir, exist_ok=True)
    overlap_file_name = output_dir / AMR_DIR_NAME / AMR_OVERLAP_FILE

    unique_amr_files = list()

    with overlap_file_name.open('w') as f:
        for i, seq in enumerate(unique_amr_seqs):
            amr_name = amr_name_from_comment(unique_amr_infos[i]["name"])
            restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
            amr_file = create_fasta_file(
                seq,
                str(amr_dir.absolute()),
                f'>{unique_amr_infos[i]["name"]}',
                restricted_amr_name
            )
            unique_amr_files.append(amr_file)
            f.write(amr_name + ":")
            if unique_amr_infos[i]["overlap_list"]:
                f.write(", ".join(e for e in unique_amr_infos[i]["overlap_list"]))
                f.write("\n")
            else:
                f.write("\n")
    return unique_amr_files
