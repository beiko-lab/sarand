"""ORF calling (pyrodigal) and annotation comparison/de-duplication primitives.

pyrodigal performs gene (CDS) calling only; it does not assign gene names or
functional products, so 'gene' and 'product' are left empty. The target AMR
gene within a neighbourhood is identified from its position (the lower-case
region of the extracted sequence), not from any annotation label.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pyrodigal

from sarand.util.logger import LOG
from sarand.util.sequence import compare_two_sequences

GeneInfo = Dict[str, Any]

# pyrodigal gene finder reused across calls; metagenomic mode (no per-sequence
# training) is required as the extracted neighbourhoods are short and varied.
_ORF_FINDER = pyrodigal.GeneFinder(meta=True)


def call_orfs(seq: str) -> List[GeneInfo]:
    """
    Call open reading frames (ORFs) in an extracted neighbourhood sequence using
    pyrodigal.

    Parameters:
        seq: the sequence to be annotated. The AMR gene is in lower case and the
            flanking neighbourhood in upper case; a trailing newline is tolerated.
    Return:
        the list of called ORFs and their details
    """
    nt_seq = seq.rstrip("\n")
    genes = _ORF_FINDER.find_genes(nt_seq.upper().encode())

    seq_info = []
    for gene in genes:
        # pyrodigal coordinates are 1-based inclusive; flip start/end on the
        # reverse strand to match the original (Bakta) output convention
        if gene.strand == 1:
            start_pos, end_pos = gene.begin, gene.end
        else:
            start_pos, end_pos = gene.end, gene.begin
        seq_info.append({
            "gene": "",
            "product": "",
            "length": (gene.end - gene.begin) + 1,
            "start_pos": start_pos,
            "end_pos": end_pos,
            "coverage": None,
            "seq_value": nt_seq,
            "seq_name": None,
            "target_amr": None,
        })
    return seq_info


def partition_genes_around_amr(
        sequence: str, seq_info: List[GeneInfo]
) -> Tuple[bool, GeneInfo | list, List[GeneInfo], List[GeneInfo], List[GeneInfo]]:
    """
    Locate the target AMR gene among the called ORFs and partition the rest into
    upstream and downstream genes.

    The AMR is the lower-case region of ``sequence``; the ORF that overlaps it
    most is flagged with ``target_amr = "yes"``.

    Return:
        (found, amr_gene, upstream_genes, downstream_genes, seq_info)
    """
    amr_start = -1
    amr_end = -1
    index = 0
    up_info = []
    down_info = []
    while index < len(sequence):
        if sequence[index].islower() and amr_start == -1:
            amr_start = index
        elif sequence[index].isupper() and amr_start > -1:
            amr_end = index - 1
            break
        index += 1
    # if there is no downstream
    if amr_end == -1 and sequence[-1].islower():
        amr_end = len(sequence) - 1
    elif amr_end == -1 or amr_start == -1:
        LOG.error("No AMR sequence (lower case string) was found in " + sequence)
        sys.exit(1)
    # find the gene that has the most overlap with the found range
    overlap_thr = 50
    found = False
    amr_info = []
    for gene_info in seq_info:
        start, end = min(gene_info["start_pos"], gene_info["end_pos"]), max(
            gene_info["start_pos"], gene_info["end_pos"]
        )
        if end < amr_start:
            up_info.append(gene_info)
        elif start > amr_end:
            down_info.append(gene_info)
        else:
            # added by 1 because in string indices start from 0
            diff = max((amr_start + 1 - start), 0) + max((end - (amr_end + 1)), 0)
            if ((1 - (float(diff) / (end - start))) * 100) > overlap_thr:
                found = True
                gene_info["target_amr"] = "yes"
                amr_info = gene_info
            elif start < amr_start:
                up_info.append(gene_info)
            else:
                down_info.append(gene_info)

    return found, amr_info, up_info, down_info, seq_info


def unnamed_genes_similar(gene_info1: GeneInfo, gene_info2: GeneInfo,
                          output_dir: str | Path, threshold: int = 90) -> bool:
    """Whether two unnamed (ORF-only) genes are significantly similar by blastn."""
    if gene_info1["gene"] != "" or gene_info2["gene"] != "":
        return False
    start1, end1 = min(gene_info1["start_pos"], gene_info1["end_pos"]), max(
        gene_info1["start_pos"], gene_info1["end_pos"]
    )
    seq1 = gene_info1["seq_value"][start1 - 1: end1 - 1]
    start2, end2 = min(gene_info2["start_pos"], gene_info2["end_pos"]), max(
        gene_info2["start_pos"], gene_info2["end_pos"]
    )
    seq2 = gene_info2["seq_value"][start2 - 1: end2 - 1]
    return compare_two_sequences(seq1, seq2, output_dir, threshold)


def annotations_identical(seq_info1: List[GeneInfo], seq_info2: List[GeneInfo],
                          out_dir: str | Path, threshold: int = 90) -> bool:
    """Whether two annotated sequences have identical gene content."""
    if len(seq_info1) == len(seq_info2):
        identical_rows = 0
        for i, gene_info1 in enumerate(seq_info1):
            gene_info2 = seq_info2[i]
            if (
                    gene_info1["gene"] == gene_info2["gene"] and gene_info1["gene"] != ""
            ) or (
                    gene_info1["gene"] == gene_info2["gene"]
                    and unnamed_genes_similar(gene_info1, gene_info2, out_dir, threshold)
            ):
                identical_rows += 1
        if identical_rows == len(seq_info1):
            return True
    return False


def similar_annotation_exists(seq_info_list: List[GeneInfo],
                              all_seq_info_lists: List[List[GeneInfo]],
                              out_dir: str | Path, threshold: int = 90) -> bool:
    """
    Whether the annotation of a new sequence already exists in the list of
    annotations extracted from other sequences (identical in 'gene', 'length',
    'start_pos' and 'end_pos' for every item).
    """
    for seq_list in all_seq_info_lists:
        if annotations_identical(seq_info_list, seq_list, out_dir, threshold):
            return True
    return False


def annotation_already_exists(seq_info_list: List[GeneInfo],
                              all_seq_info_lists: List[List[GeneInfo]],
                              out_dir: str | Path) -> bool:
    """Like ``similar_annotation_exists`` but requiring an exact (100%) match."""
    for seq_list in all_seq_info_lists:
        if annotations_identical(seq_info_list, seq_list, out_dir, 100):
            return True
    return False
