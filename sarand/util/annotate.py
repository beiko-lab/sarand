"""ORF calling (pyrodigal) and annotation comparison/de-duplication primitives.

pyrodigal performs gene (CDS) calling only; it does not assign gene names or
functional products, so 'gene' and 'product' are left empty. The target AMR
gene within a neighborhood is identified from its position (the lower-case
region of the extracted sequence), not from any annotation label.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pyrodigal

from sarand.util.logger import LOG
from sarand.util.sequence import similar_sequence_pairs

GeneInfo = Dict[str, Any]

# pyrodigal gene finder reused across calls; metagenomic mode (no per-sequence
# training) is required as the extracted neighborhoods are short and varied.
_ORF_FINDER = pyrodigal.GeneFinder(meta=True)


def call_orfs(seq: str) -> List[GeneInfo]:
    """
    Call open reading frames (ORFs) in an extracted neighborhood sequence using
    pyrodigal.

    Parameters:
        seq: the sequence to be annotated. The AMR gene is in lower case and the
            flanking neighborhood in upper case; a trailing newline is tolerated.
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
            # retain pyrodigal's ORF output so the called ORFs can be written
            # out later (e.g. for the coverage-filtered neighborhoods)
            "strand": gene.strand,
            "partial": f"{int(gene.partial_begin)}{int(gene.partial_end)}",
            "nt_seq": gene.sequence(),
            "aa_seq": gene.translate(),
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


def _orf_nt_seq(gene_info: GeneInfo) -> str:
    """The nucleotide sub-sequence of a called ORF within its neighborhood."""
    start = min(gene_info["start_pos"], gene_info["end_pos"])
    end = max(gene_info["start_pos"], gene_info["end_pos"])
    return gene_info["seq_value"][start - 1: end - 1]


def similar_annotation_exists(seq_info_list: List[GeneInfo],
                              all_seq_info_lists: List[List[GeneInfo]],
                              out_dir: str | Path, threshold: int = 90) -> bool:
    """
    Whether a neighborhood's ORF annotation already exists among previously kept
    neighborhoods (i.e. it is a near-duplicate to collapse).

    Two neighborhoods are considered identical when they have the same number of
    ORFs and, at every position, the ORFs either share a non-empty gene name or
    (for unnamed ORFs) have near-identical nucleotide sequences. Every ORF
    sequence comparison the check might need is performed in a single batched
    minimap2 alignment rather than one alignment per pair.
    """
    # only neighborhoods with the same number of ORFs can be identical
    candidates = [
        (idx, kept) for idx, kept in enumerate(all_seq_info_lists)
        if len(kept) == len(seq_info_list)
    ]
    if not candidates:
        return False

    # collect every unnamed-ORF comparison the check might need, keyed so each
    # (query, subject) result maps back to a (position, candidate) pair, and run
    # them all in one minimap2 alignment
    queries = {
        f"q{i}": _orf_nt_seq(g)
        for i, g in enumerate(seq_info_list) if g["gene"] == ""
    }
    subjects = {
        f"s{idx}_{i}": _orf_nt_seq(g)
        for idx, kept in candidates
        for i, g in enumerate(kept) if g["gene"] == ""
    }
    matched = similar_sequence_pairs(queries, subjects, threshold)

    for idx, kept in candidates:
        for i, gene_info in enumerate(seq_info_list):
            other = kept[i]
            if gene_info["gene"] != "" or other["gene"] != "":
                # named ORF(s): identical only if they share a non-empty name
                if not (gene_info["gene"] == other["gene"] and gene_info["gene"] != ""):
                    break
            elif (f"q{i}", f"s{idx}_{i}") not in matched:
                # both unnamed: require near-identical nucleotide sequences
                break
        else:
            # every position matched -> this neighborhood is a duplicate
            return True
    return False


def write_orf_files(seq_info_lists: List[List[GeneInfo]],
                    out_prefix: str | Path) -> Tuple[Path, Path, Path]:
    """Write the pyrodigal ORFs of each neighborhood to FASTA and GFF3 files.

    Parameters:
        seq_info_lists: one entry per neighborhood, each a list of ORF dicts as
            produced by ``call_orfs`` (so each carries ``nt_seq``/``aa_seq``).
        out_prefix: path prefix; ``<prefix>.ffn`` (nucleotide), ``<prefix>.faa``
            (protein) and ``<prefix>.gff`` (GFF3) are written.
    Return:
        the (ffn, faa, gff) paths written.
    """
    out_prefix = Path(out_prefix)
    ffn_path = out_prefix.with_suffix(".ffn")
    faa_path = out_prefix.with_suffix(".faa")
    gff_path = out_prefix.with_suffix(".gff")
    with open(ffn_path, "w") as ffn, open(faa_path, "w") as faa, \
            open(gff_path, "w") as gff:
        gff.write("##gff-version 3\n")
        for seq_info in seq_info_lists:
            for orf_index, gene_info in enumerate(seq_info, start=1):
                seqid = str(gene_info["seq_name"])
                orf_id = f"{seqid}_orf{orf_index}"
                lo = min(gene_info["start_pos"], gene_info["end_pos"])
                hi = max(gene_info["start_pos"], gene_info["end_pos"])
                strand = "+" if gene_info.get("strand", 1) == 1 else "-"
                partial = gene_info.get("partial", "00")
                is_target = gene_info.get("target_amr") == "yes"
                header = f"{orf_id} # {lo} # {hi} # {strand} # partial={partial}"
                ffn.write(f">{header}\n{gene_info.get('nt_seq', '')}\n")
                faa.write(f">{header}\n{gene_info.get('aa_seq', '')}\n")
                gff.write(
                    f"{seqid}\tpyrodigal\tCDS\t{lo}\t{hi}\t.\t{strand}\t0\t"
                    f"ID={orf_id};partial={partial};target={'yes' if is_target else 'no'}\n"
                )
    return ffn_path, faa_path, gff_path
