"""Stage 3 of the graph pipeline: ORF-annotate the extracted neighbourhoods.

Each extracted neighbourhood sequence is run through pyrodigal ORF calling, the
target AMR is located within it, per-gene coverage is computed, and the results
are written to the annotation CSVs.
"""
import csv
import shutil
import sys
from multiprocessing.pool import Pool
from pathlib import Path

from sarand.config import ANNOTATION_DIR, SEQ_NAME_PREFIX
from sarand.coverage import find_gene_coverage, read_path_coverage_info
from sarand.util.annotate import (
    annotation_already_exists,
    call_orfs,
    partition_genes_around_amr,
)
from sarand.util.file import try_dump_to_disk
from sarand.util.logger import LOG
from sarand.util.naming import extract_name_from_file_name
from sarand.util.sequence import retrieve_AMR


def write_annotation_row(annotation_writer, visual_annotation_writer, gene_info, found,
                         len_seq=None):
    """
    Write one annotation row to the detail file (and the unique/trimmed file when
    the annotation has not been seen before).
    Parameters:
        annotation_writer:	annotation file containing all annotations
        visual_annotation_writer: annotation file containing unique annotations
        gene_info: annotation info
        found: if True, the annotation info has already been found in other sequences
    """
    seq = gene_info["seq_value"]
    seq_description = gene_info["seq_name"]
    if len_seq is None:
        len_seq = len(seq)
    row = [
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
    ]
    annotation_writer.writerow(row)
    if not found:
        visual_annotation_writer.writerow(row)


def annotate_one_sequence(seq_pair):
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
        amr_name,
        path_info_file,
        neighborhood_seq_file,
        annotate_dir,
        core_num,
        annotation_writer,
        trimmed_annotation_writer,
        gene_file,
        error_file,
):
    """
    Annotate (in parallel) the neighbourhood sequences extracted for one AMR.
    Parameters:
        amr_name: the name of AMR
        path_info_file: the information of nodes representing the extracted sequence
        neighborhood_seq_file: the file containing extracted sequences
        annotate_dir: the directory to store annotation results
        core_num: the number of cores for parallel processing
        annotation_writer: the file to store annotation results
        trimmed_annotation_writer: the file to store unique annotation results
        gene_file: the file to store gene names in annotation
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
    LOG.debug("Reading " + neighborhood_seq_file + " for " + amr_name)
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
        # extract amr from seq_info
        amr_found, amr_info, up_info, down_info, seq_info = partition_genes_around_amr(
            line[:-1], seq_info
        )
        if not amr_found:
            LOG.error("ERROR: no target amr was found in the extracted sequence")
            error_writer.write(
                amr_name
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
        # Check if this annotation has already been found
        found = annotation_already_exists(seq_info, all_seq_info_list, annotate_dir)
        if not found:
            all_seq_info_list.append(seq_info)
        myLine1 = seq_description + ":\t"
        # write annotation info into the files
        for j, gene_info in enumerate(seq_info):
            coverage = coverage_list[j] if coverage_list else -1
            gene_info["coverage"] = coverage
            gene_info["seq_name"] = seq_description
            write_annotation_row(
                annotation_writer, trimmed_annotation_writer, gene_info, found
            )
            if gene_info["gene"] == "":
                myLine1 += "UNKNOWN---"
            else:
                myLine1 += gene_info["gene"] + "---"
        gene_file.write(myLine1[:-3] + "\n")
    if not all_seq_info_list:
        error_writer.write(amr_name + " no annotation was found in the graph.\n")
    error_writer.close()
    return all_seq_info_list


def annotate_neighborhood(
        amr_name,
        neighborhood_seq_file,
        path_info_file,
        seq_length,
        output_dir,
        output_name="",
        core_num=4,
):
    """
    Annotate the neighbourhood sequences extracted from the assembly graph for one
    AMR and summarise the results into the annotation CSV files.
    Parameters:
        amr_name:	the name of target AMR
        neighborhood_seq_file:	the file containing all extracted neighbourhood sequences
        path_info_file: the per-node path/coverage info file
        seq_length:	the length of neighbourhood sequence extracted from each side
        output_dir:	the path for the output directory
        output_name: the suffix used to distinguish output files (usually the AMR name)
        core_num: the number of cores for parallel processing
    Return:
        (all_seq_info_list, trimmed_annotation_info_name)
    """
    LOG.info("Annotating " + amr_name)
    # initializing required files and directories
    annotations_dir = (
        Path(output_dir) / ANNOTATION_DIR / (ANNOTATION_DIR + "_" + str(seq_length))
    )
    annotate_dir = annotations_dir / ("annotation" + output_name + "_" + str(seq_length))
    if annotate_dir.exists():
        try:
            shutil.rmtree(annotate_dir)
        except OSError as e:
            LOG.error("Error: %s - %s." % (e.filename, e.strerror))
    annotate_dir.mkdir(parents=True)
    error_file = annotations_dir / "not_found_annotation_amrs_in_graph.txt"
    annotation_detail_name = annotate_dir / ("annotation_detail" + output_name + ".csv")
    trimmed_annotation_info_name = annotate_dir / (
        "trimmed_annotation_info" + output_name + ".csv"
    )
    annotation_detail = open(annotation_detail_name, mode="w", newline="")
    trimmed_annotation_info = open(trimmed_annotation_info_name, mode="w", newline="")
    annotation_writer = csv.writer(annotation_detail)
    trimmed_annotation_writer = csv.writer(trimmed_annotation_info)
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
    write_annotation_row(
        annotation_writer,
        trimmed_annotation_writer,
        header,
        False,
        "seq_length",
    )
    gene_file_name = annotate_dir / ("seq_comparison_genes" + output_name + ".txt")
    gene_file = open(gene_file_name, "w")

    # annotate the sequences extracted from the assembly graph
    all_seq_info_list = annotate_graph_sequences(
        amr_name,
        path_info_file,
        neighborhood_seq_file,
        annotate_dir,
        core_num,
        annotation_writer,
        trimmed_annotation_writer,
        gene_file,
        error_file,
    )
    LOG.debug(
        f"The comparison of neighborhood sequences are available in "
        f"{annotation_detail_name}, {gene_file_name}"
    )
    annotation_detail.close()
    trimmed_annotation_info.close()
    gene_file.close()

    return all_seq_info_list, str(trimmed_annotation_info_name)


def find_seq_and_path_files(amr_name, sequences_file_names, path_info_file_names,
                            seq_length):
    """
    Return the extracted-sequence file and the path-info file produced for a given
    AMR at a given neighbourhood length (or -1 if not found).
    """
    seq_file = -1
    for file_name in sequences_file_names:
        if SEQ_NAME_PREFIX + amr_name + "_" + str(seq_length) in file_name:
            seq_file = file_name
    path_file = -1
    for file_name in path_info_file_names:
        if SEQ_NAME_PREFIX + amr_name + "_" + str(seq_length) in file_name:
            path_file = file_name
    return seq_file, path_file


def annotate_all_amrs(params, seq_files, path_info_files, amr_files, debug: bool):
    """
    Annotate the neighbourhood sequences of every AMR.
    Parameters:
        params: the parsed CLI parameters
        seq_files: the list of neighbourhood sequence files
        path_info_files: the list of per-node path/coverage files
        amr_files: the list of files containing AMRs
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
    for amr_file in amr_files:
        restricted_amr_name = extract_name_from_file_name(amr_file)
        _, amr_name = retrieve_AMR(amr_file)
        neighborhood_file, nodes_info_file = find_seq_and_path_files(
            restricted_amr_name,
            neighborhood_files,
            nodes_info_files,
            params.neighbourhood_length,
        )
        if neighborhood_file == -1:
            LOG.error(
                "no sequence file for "
                + amr_file
                + " was found! We looked for a file like "
                + restricted_amr_name
            )
            sys.exit(1)
        all_seq_info_list, annotation_file = annotate_neighborhood(
            amr_name,
            neighborhood_file,
            nodes_info_file,
            params.neighbourhood_length,
            params.output_dir,
            "_" + restricted_amr_name,
            params.num_cores,
        )
        all_seq_info_lists.append(all_seq_info_list)
        annotation_files.append(annotation_file)

    if debug:
        try_dump_to_disk(
            {'all_seq_info_lists': all_seq_info_lists, 'annotation_files': annotation_files},
            Path(params.output_dir) / 'sequences_info' / 'debug_seq_annotation_main.json'
        )

    return all_seq_info_lists, annotation_files
