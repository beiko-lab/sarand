"""Per-gene coverage calculation and coverage-consistency filtering of annotations."""
from __future__ import annotations

import argparse
import copy
import csv
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

from sarand.config import ANNOTATION_DIR
from sarand.util.annotate import similar_annotation_exists, write_orf_files
from sarand.util.logger import LOG
from sarand.util.naming import extract_name_from_file_name

GeneInfo = Dict[str, Any]
NodeInfo = Dict[str, Any]


def read_path_coverage_info(path_info_file: str | Path) -> List[List[NodeInfo]]:
    """
    Read the csv file listing the nodes (and their coverage) that make up each
    extracted sequence.
    Parameter:
        path_info_file: the csv file
    Return:
        stored info in a list
    """
    seq_info_list = []
    seq_info = []
    with open(path_info_file, "r") as myfile:
        myreader = csv.DictReader(myfile)
        old_seq = ""
        for row in myreader:
            # some path info files repeat the header / duplicate rows; skip them
            if row["coverage"] == "coverage":
                continue
            node_info = {
                "node": row["node"],
                "coverage": float(row["coverage"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
            }
            cur_seq = row["sequence"]
            if cur_seq != old_seq:
                if seq_info:
                    seq_info_list.append(seq_info)
                seq_info = []
                old_seq = cur_seq
            seq_info.append(node_info)
        seq_info_list.append(seq_info)
    return seq_info_list


def find_gene_coverage(seq_info_list: List[GeneInfo],
                       path_info: List[NodeInfo]) -> List[float]:
    """
    Calculate the coverage of genes available in a sequence based on the coverage
    of the nodes representing them
    Parameter:
        seq_info_list: the list of genes (and their info) annotated in a sequence
        path_info:	the list of nodes (and their info) representing the sequence
    Return:
        the list of calculated gene-coverages
    """
    coverage_list = []
    for seq_info in seq_info_list:
        sum_coverage = 0
        # minus 1 because the ORF caller returns 1-based start positions
        start, end = seq_info["start_pos"] - 1, seq_info["end_pos"] - 1
        if start > end:
            start, end = end, start
        found = False
        for j, node_info in enumerate(path_info):
            if node_info["start"] <= start and node_info["end"] >= start:
                found = True
                # the annotated gene is represented by a single node
                if node_info["end"] >= end:
                    sum_coverage = (end - start + 1) * node_info["coverage"]
                else:
                    # calculate the length of node representing the gene
                    sum_coverage = (node_info["end"] - start + 1) * node_info[
                        "coverage"
                    ]
                    n_index = j + 1
                    assert n_index < len(
                        path_info
                    ), "wrong index calculated for path_info!!"
                    while (
                            n_index < len(path_info) and path_info[n_index]["start"] <= end
                    ):
                        if path_info[n_index]["end"] >= end:
                            sum_coverage += (
                                                    end - path_info[n_index]["start"] + 1
                                            ) * path_info[n_index]["coverage"]
                            break
                        else:
                            sum_coverage += (
                                                    path_info[n_index]["end"]
                                                    - path_info[n_index]["start"]
                                                    + 1
                                            ) * path_info[n_index]["coverage"]
                        n_index += 1
                coverage_list.append(sum_coverage / (end - start + 1))
                break
        if not found:
            LOG.error("ERROR: no nodes were found for this AMR gene!!!")
    return coverage_list


def find_amr_coverage_in_seq(seq_info: List[GeneInfo]) -> Tuple[float, int, bool]:
    """
    Find the target AMR within an annotated sequence and its coverage.

    The AMR is the lower-case region of the extracted sequence; the annotated
    gene that overlaps it most is taken as the AMR.
    Return:
        the coverage of the AMR, its index in seq_info, and False if a gene was
        found in the lower-case area of the sequence; otherwise True.
    """
    error = True
    amr_coverage = 0
    # find the indeces of lower case string (amr sequence) in the extracted sequence
    sequence = seq_info[0]["seq_value"]
    amr_start = -1
    amr_end = -1
    index = 0
    while index < len(sequence):
        if sequence[index].islower() and amr_start == -1:
            amr_start = index
        elif sequence[index].isupper() and amr_start > -1:
            amr_end = index - 1
            break
        index += 1
    # find the gene with the most overlap with the found range, above most_overlap
    most_overlap = 50
    amr_index = -1
    for i, gene_info in enumerate(seq_info):
        start, end = min(gene_info["start_pos"], gene_info["end_pos"]), max(
            gene_info["start_pos"], gene_info["end_pos"]
        )
        if end < amr_start or start > amr_end:
            continue
        else:
            # added by 1 because in string indices start from 0
            diff = max((amr_start + 1 - start), 0) + max((end - (amr_end + 1)), 0)
            if ((1 - (float(diff) / (end - start))) * 100) > most_overlap:
                most_overlap = (1 - (float(diff) / (end - start))) * 100
                amr_coverage = gene_info["coverage"]
                amr_index = i
                error = False
    return amr_coverage, amr_index, error


def filter_by_coverage_consistency(
        seq_info_list_input: List[List[GeneInfo]], coverage_thr: int,
        amr_name: str, annotate_dir: str | Path,
) -> Tuple[str | Path, int]:
    """
    Compare the coverage of each gene in the neighborhood with that of the AMR
    and drop genes (plus everything upstream/downstream of them) whose coverage
    differs from the AMR's by more than ``coverage_thr``. Writes the surviving
    annotations to ``coverage_annotation_<thr>_<amr_name>.csv``.
    """
    seq_info_list = copy.deepcopy(seq_info_list_input)
    # extract amr info
    amr_coverages = []
    amr_indeces = []
    for seq_info in seq_info_list:
        found_amr = False
        for gene_counter, gene_info in enumerate(seq_info):
            if gene_info["coverage"] is None:
                LOG.info("Coverage information are not available for " + amr_name)
                return "", []
            coverage = round(gene_info["coverage"], 2)
            if gene_info["target_amr"] == "yes":
                amr_coverages.append(coverage)
                amr_indeces.append(gene_counter)
                found_amr = True
                break
        if not found_amr:
            (
                amr_coverage,
                amr_index,
                error,
            ) = find_amr_coverage_in_seq(seq_info)
            if error:
                LOG.error(
                    "No target amr was found for "
                    + str(seq_info)
                    + " regarding "
                    + amr_name
                )
                sys.exit(1)
            else:
                amr_coverages.append(amr_coverage)
                amr_indeces.append(amr_index)
    if len(amr_coverages) != len(seq_info_list):
        LOG.error(
            "Inconsistency between the number of sequences and found amr-coverages!"
        )
        sys.exit(1)
    # remove genes with inconsistent coverage and whatever comes before them if
    # upstream OR after them if downstream
    remained_seqs = []
    for i, seq_info in enumerate(seq_info_list):
        # find the genes need to be removed
        to_be_removed_genes = []
        for j, gene_info in enumerate(seq_info):
            if abs(gene_info["coverage"] - amr_coverages[i]) > coverage_thr:
                if j < amr_indeces[i]:
                    for k in range(j + 1):
                        if k not in to_be_removed_genes:
                            to_be_removed_genes.append(k)
                elif j > amr_indeces[i]:
                    for k in range(j, len(seq_info)):
                        if k not in to_be_removed_genes:
                            to_be_removed_genes.append(k)
                    break
        for j in reversed(range(len(seq_info))):
            if j in to_be_removed_genes:
                del seq_info[j]
        # check if the remained sequence already exists in the seq_info_list
        if seq_info and not similar_annotation_exists(
                seq_info, remained_seqs, annotate_dir
        ):
            remained_seqs.append(seq_info)

    # Initialize coverage file
    coverage_annotation = Path(annotate_dir) / (
        "coverage_annotation_" + str(coverage_thr) + "_" + amr_name + ".csv"
    )
    with open(coverage_annotation, "w") as fd:
        writer = csv.writer(fd)
        writer.writerow(
            [
                "seq_name",
                "seq_value",
                "seq_length",
                "coverage",
                "length",
                "start_pos",
                "end_pos",
                "target_amr",
            ]
        )
        # write the ORFs of extracted sequences with consistent coverage
        for seq_info in remained_seqs:
            for gene_info in seq_info:
                writer.writerow(
                    [
                        gene_info["seq_name"],
                        gene_info["seq_value"],
                        len(gene_info["seq_value"]),
                        gene_info["coverage"],
                        gene_info["length"],
                        gene_info["start_pos"],
                        gene_info["end_pos"],
                        gene_info["target_amr"],
                    ]
                )

    return coverage_annotation, remained_seqs


def _gene_path(seq_info: List[GeneInfo]) -> str:
    """Return the ordered ORF coordinates of a neighborhood as 'lo-hi;lo-hi;...'."""
    return ";".join(
        f"{min(g['start_pos'], g['end_pos'])}-{max(g['start_pos'], g['end_pos'])}"
        for g in seq_info
    )


def _format_coverage(coverage: Any) -> str:
    """Render a coverage value for the combined CSV (rounded; blank if missing)."""
    if coverage is None or coverage == "":
        return ""
    return str(round(float(coverage), 2))


def _target_coverage(seq_info: List[GeneInfo]) -> Any:
    """Return the coverage of the target ORF (the one flagged target_amr='yes')."""
    for gene_info in seq_info:
        if gene_info.get("target_amr") == "yes":
            return gene_info["coverage"]
    return ""


def _target_gene_coords(seq_info: List[GeneInfo]) -> str:
    """Return the 'lo-hi' coordinates of the target ORF (flagged target_amr='yes')."""
    for gene_info in seq_info:
        if gene_info.get("target_amr") == "yes":
            return (
                f"{min(gene_info['start_pos'], gene_info['end_pos'])}"
                f"-{max(gene_info['start_pos'], gene_info['end_pos'])}"
            )
    return ""


def write_combined_final_neighborhoods(
        final_neighborhoods: List[Tuple[str, List[GeneInfo]]],
        fasta_file: str | Path,
        csv_file: str | Path,
) -> None:
    """
    Write the combined outputs for the final (coverage-filtered) neighborhoods:
    a single FASTA of the neighborhood sequences and a single CSV summarising the
    target gene, ORF path and coverages of each.
    Parameters:
        final_neighborhoods: (target_name, seq_info) pairs, one per neighborhood
        fasta_file: path of the combined FASTA to write
        csv_file: path of the combined summary CSV to write
    """
    with open(fasta_file, "w") as fd:
        for target_name, seq_info in final_neighborhoods:
            seq_name = seq_info[0]["seq_name"]
            fd.write(f">{target_name}_{seq_name}\n{seq_info[0]['seq_value']}\n")

    with open(csv_file, "w", newline="") as fd:
        writer = csv.writer(fd)
        writer.writerow(
            ["target_name", "seq_name", "target_gene", "gene_path",
             "target_coverage", "coverages"]
        )
        for target_name, seq_info in final_neighborhoods:
            writer.writerow([
                target_name,
                seq_info[0]["seq_name"],
                _target_gene_coords(seq_info),
                _gene_path(seq_info),
                _format_coverage(_target_coverage(seq_info)),
                ";".join(_format_coverage(g["coverage"]) for g in seq_info),
            ])


def trim_annotations_by_coverage(params: argparse.Namespace, target_files: List[str],
                                 all_seq_info_lists: List[List[List[GeneInfo]]]) -> List[str | Path]:
    """
    Determine the final neighborhoods for each target and write the combined
    outputs under ``final_neighborhoods``.

    When ``--coverage_difference`` is positive each target's neighborhood
    annotations are filtered by coverage consistency (writing a per-target
    coverage_annotation csv); otherwise every annotated neighborhood is kept as
    final. For each target the called ORFs of its final neighborhoods are written
    (``orfs_<target>.{ffn,faa,gff}``), and across all targets a single combined
    FASTA of the final neighborhood sequences and a single combined summary CSV
    are written.
    Parameters:
        params: the parsed CLI parameters
        target_files: the list of files containing target genes
        all_seq_info_lists: the annotations of the neighborhood sequences
    """
    annotations_dir = Path(params.output_dir) / ANNOTATION_DIR
    coverage_annotation_list: List[str | Path] = []
    final_neighborhoods: List[Tuple[str, List[GeneInfo]]] = []
    for i, target_file in enumerate(target_files):
        restricted_target_name = extract_name_from_file_name(target_file)
        annotate_dir = annotations_dir / (
            "annotation_" + restricted_target_name + "_" + str(params.neighborhood_length)
        )
        if params.coverage_difference > 0:
            # remove extracted sequences with inconsistent coverage
            coverage_annotation, remained_seqs = filter_by_coverage_consistency(
                all_seq_info_lists[i],
                params.coverage_difference,
                restricted_target_name,
                annotate_dir,
            )
            coverage_annotation_list.append(coverage_annotation)
        else:
            # no coverage threshold: every annotated neighborhood is final
            remained_seqs = all_seq_info_lists[i]

        # write the pyrodigal ORFs called per final neighborhood
        write_orf_files(remained_seqs, annotate_dir / ("orfs_" + restricted_target_name))

        for seq_info in remained_seqs:
            if seq_info:
                final_neighborhoods.append((restricted_target_name, seq_info))

    write_combined_final_neighborhoods(
        final_neighborhoods,
        annotations_dir / "final_neighborhoods.fasta",
        annotations_dir / "final_neighborhoods.csv",
    )
    return coverage_annotation_list
