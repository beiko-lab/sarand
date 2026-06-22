"""Stage 3 of the graph pipeline: ORF-annotate the extracted neighborhoods.

Each extracted neighborhood sequence is run through pyrodigal ORF calling, the
target gene is located within it, per-gene coverage is computed, and the results
are written to the annotation CSV.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from sarand.config import ANNOTATION_DIR, SEQ_NAME_PREFIX

GeneInfo = Dict[str, Any]
from sarand.coverage import find_gene_coverage, read_path_coverage_info
from sarand.util.annotate import (
    call_orfs,
    partition_genes_around_amr,
)
from sarand.util.file import try_dump_to_disk
from sarand.util.logger import LOG
from sarand.util.naming import extract_name_from_file_name
from sarand.util.sequence import retrieve_target


def write_annotation_row(annotation_writer: Any, gene_info: GeneInfo,
                         len_seq: Optional[int] = None) -> None:
    """
    Write one annotation row to the detail file.
    Parameters:
        annotation_writer: annotation file containing all annotations
        gene_info: annotation info
    """
    seq = gene_info["seq_value"]
    seq_description = gene_info["seq_name"]
    if len_seq is None:
        len_seq = len(seq)
    annotation_writer.writerow([
        seq_description,
        seq,
        len_seq,
        gene_info["gene"],
        gene_info["product"],
        gene_info["length"],
        gene_info["start_pos"],
        gene_info["end_pos"],
        gene_info["coverage"],
        gene_info["target_amr"],
    ])


def annotate_one_sequence(seq_pair: Tuple[int, str]) -> List[GeneInfo]:
    """
    Worker for parallel annotation: call ORFs in a single extracted sequence.
    Parameters:
        seq_pair: the (index, sequence) pair to be annotated
    Return:
        the list of called ORFs
    """
    counter, ext_seq = seq_pair
    return call_orfs(ext_seq)


def annotate_graph_sequences(
        target_name: str,
        path_info_file: str | Path | int,
        neighborhood_seq_file: str | Path,
        core_num: int,
        annotation_writer: Any,
        error_file: str | Path,
) -> List[List[GeneInfo]]:
    """
    Annotate (in parallel) the neighborhood sequences extracted for one target.
    Parameters:
        target_name: the name of the target gene
        path_info_file: the information of nodes representing the extracted sequence
        neighborhood_seq_file: the file containing extracted sequences
        core_num: the number of cores for parallel processing
        annotation_writer: the file to store annotation results
        error_file: the file to store errors
    Return:
        the list of annotated genes and their details
    """
    error_writer = open(error_file, "a")
    # Read path_info from the file
    path_info_list = []
    if path_info_file != -1:
        path_info_list = read_path_coverage_info(path_info_file)
    # find the list of all extracted sequences
    LOG.debug("Reading " + neighborhood_seq_file + " for " + target_name)
    sequence_list = []
    counter = 1
    with open(neighborhood_seq_file, "r") as read_obj:
        for line in read_obj:
            if (
                    line.startswith(">")
                    or line.startswith("Path")
                    or line.startswith("The")
            ):
                continue
            sequence_list.append((counter, line))
            counter += 1

    # Parallel annotation (skip the pool overhead when single-threaded)
    if core_num == 1:
        seq_info_list = list()
        for x in sequence_list:
            seq_info_list.append(annotate_one_sequence(x))
    else:
        with Pool(core_num) as p:
            seq_info_list = p.map(annotate_one_sequence, sequence_list)

    # Further processing of result of parallel annotation
    all_seq_info_list = []
    for i, seq_pair in enumerate(sequence_list):
        counter, line = seq_pair
        seq_description = "extracted" + str(counter)
        seq_info = seq_info_list[i]
        if not seq_info or len(seq_info) == 0:
            continue
        # locate the target gene among the called ORFs
        target_found, target_info, up_info, down_info, seq_info = partition_genes_around_amr(
            line[:-1], seq_info
        )
        if not target_found:
            LOG.error("ERROR: no target gene was found in the extracted sequence")
            error_writer.write(
                target_name
                + " annotation not found! "
                + " seq_info: "
                + str(seq_info)
                + "\n"
            )
            continue
        # calculate coverage for the genes available in the annotation
        coverage_list = []
        if path_info_list:
            coverage_list = find_gene_coverage(seq_info, path_info_list[counter - 1])
        all_seq_info_list.append(seq_info)
        # write annotation info into the file
        for j, gene_info in enumerate(seq_info):
            coverage = coverage_list[j] if coverage_list else -1
            gene_info["coverage"] = coverage
            gene_info["seq_name"] = seq_description
            write_annotation_row(annotation_writer, gene_info)
    if not all_seq_info_list:
        error_writer.write(target_name + " no annotation was found in the graph.\n")
    error_writer.close()
    return all_seq_info_list


def annotate_neighborhood(
        target_name: str,
        neighborhood_seq_file: str | Path,
        path_info_file: str | Path | int,
        seq_length: int,
        output_dir: str | Path,
        output_name: str = "",
        core_num: int = 4,
) -> Tuple[List[List[GeneInfo]], str]:
    """
    Annotate the neighborhood sequences extracted from the assembly graph for one
    target and summarise the results into the annotation CSV.
    Parameters:
        target_name: the name of the target gene
        neighborhood_seq_file: the file containing all extracted neighborhood sequences
        path_info_file: the per-node path/coverage info file
        seq_length: the length of neighborhood sequence extracted from each side
        output_dir: the path for the output directory
        output_name: the suffix used to distinguish output files (usually the target name)
        core_num: the number of cores for parallel processing
    Return:
        (all_seq_info_list, annotation_detail_name)
    """
    LOG.info("Annotating " + target_name)
    # initializing required files and directories
    annotations_dir = Path(output_dir) / ANNOTATION_DIR
    annotate_dir = annotations_dir / ("annotation" + output_name + "_" + str(seq_length))
    if annotate_dir.exists():
        try:
            shutil.rmtree(annotate_dir)
        except OSError as e:
            LOG.error("Error: %s - %s." % (e.filename, e.strerror))
    annotate_dir.mkdir(parents=True)
    error_file = annotations_dir / "not_found_annotation_targets_in_graph.txt"
    annotation_detail_name = annotate_dir / ("annotation_detail" + output_name + ".csv")
    annotation_detail = open(annotation_detail_name, mode="w", newline="")
    annotation_writer = csv.writer(annotation_detail)
    header = {
        "seq_value": "seq_value",
        "gene": "gene",
        "product": "product",
        "length": "length",
        "start_pos": "start_pos",
        "end_pos": "end_pos",
        "coverage": "coverage",
        "seq_name": "seq_name",
        "target_amr": "target_amr",
    }
    write_annotation_row(annotation_writer, header, "seq_length")

    # annotate the sequences extracted from the assembly graph
    all_seq_info_list = annotate_graph_sequences(
        target_name,
        path_info_file,
        neighborhood_seq_file,
        core_num,
        annotation_writer,
        error_file,
    )
    LOG.debug(f"The neighborhood annotations are available in {annotation_detail_name}")
    annotation_detail.close()

    return all_seq_info_list, str(annotation_detail_name)


def find_seq_and_path_files(
        target_name: str,
        sequences_file_names: List[str],
        path_info_file_names: List[str],
        seq_length: int,
) -> Tuple[str | int, str | int]:
    """
    Return the extracted-sequence file and the path-info file produced for a given
    target at a given neighborhood length (or -1 if not found).
    """
    seq_file = -1
    for file_name in sequences_file_names:
        if SEQ_NAME_PREFIX + target_name + "_" + str(seq_length) in file_name:
            seq_file = file_name
    path_file = -1
    for file_name in path_info_file_names:
        if SEQ_NAME_PREFIX + target_name + "_" + str(seq_length) in file_name:
            path_file = file_name
    return seq_file, path_file


def annotate_all_neighborhoods(
        params: argparse.Namespace,
        seq_files: List[str],
        path_info_files: List[str],
        target_files: List[str],
        debug: bool,
) -> Tuple[List[List[List[GeneInfo]]], List[str]]:
    """
    Annotate the neighborhood sequences of every target.
    Parameters:
        params: the parsed CLI parameters
        seq_files: the list of neighborhood sequence files
        path_info_files: the list of per-node path/coverage files
        target_files: the list of files containing target genes
        debug: True if additional files should be created, False otherwise.
    Return:
        (all_seq_info_lists, annotation_files)
    """
    LOG.info("Neighborhood Annotation...")
    if seq_files:
        neighborhood_files = seq_files
    else:
        LOG.error("No file containing the extracted neighborhood sequences is available!")
        sys.exit(1)

    if path_info_files:
        nodes_info_files = path_info_files
    else:
        LOG.error("No file containing path info for neighborhood sequences is available!")
        sys.exit(1)

    all_seq_info_lists = []
    annotation_files = []
    for target_file in target_files:
        restricted_target_name = extract_name_from_file_name(target_file)
        _, target_name = retrieve_target(target_file)
        neighborhood_file, nodes_info_file = find_seq_and_path_files(
            restricted_target_name,
            neighborhood_files,
            nodes_info_files,
            params.neighborhood_length,
        )
        if neighborhood_file == -1:
            LOG.error(
                "no sequence file for "
                + target_file
                + " was found! We looked for a file like "
                + restricted_target_name
            )
            sys.exit(1)
        all_seq_info_list, annotation_file = annotate_neighborhood(
            target_name,
            neighborhood_file,
            nodes_info_file,
            params.neighborhood_length,
            params.output_dir,
            "_" + restricted_target_name,
            params.num_cores,
        )
        all_seq_info_lists.append(all_seq_info_list)
        annotation_files.append(annotation_file)

    if debug:
        try_dump_to_disk(
            {'all_seq_info_lists': all_seq_info_lists, 'annotation_files': annotation_files},
            Path(params.output_dir) / 'debug_seq_annotation_main.json'
        )

    return all_seq_info_lists, annotation_files
