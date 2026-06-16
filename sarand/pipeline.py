"""Top-level orchestration of the main (assembly-graph) sarand pipeline.

Stages: find target genes in the graph -> extract their neighbourhoods -> ORF
annotate the neighbourhoods -> filter the annotations by coverage consistency.
"""
from __future__ import annotations

import argparse
import sys

from sarand.target_finder import find_all_targets_in_graph
from sarand.annotation import annotate_all_neighbourhoods
from sarand.coverage import trim_annotations_by_coverage
from sarand.extract_neighborhood import extract_target_neighbourhoods
from sarand.util.logger import LOG


def run_graph_pipeline(params: argparse.Namespace) -> None:
    """Run the full assembly-graph sarand workflow."""
    LOG.info("Starting analysis...")

    # Stage 1: find the target genes in the assembly graph
    # this uses Bandage's BLAST implementation for graphs
    LOG.info(f"Finding target genes {params.target_genes} in the assembly graph: {params.input_gfa}")
    unique_target_files, unique_target_path_list = find_all_targets_in_graph(
        params.input_gfa,
        params.output_dir,
        params.target_genes,
        params.min_target_identity,
        params.min_target_coverage,
        params.num_cores,
        params.keep_intermediate_files,
        params.debug,
    )

    if not unique_target_files:
        LOG.error("No target genes were found in graph!")
        sys.exit(1)

    if not (unique_target_path_list and len(unique_target_path_list) == len(unique_target_files)):
        LOG.error("Target alignment data is missing/inconsistent")
        sys.exit(1)

    # pair each found target with its alignment info
    target_seq_align_info = []
    for i, target_file in enumerate(unique_target_files):
        target_seq_align_info.append((target_file, unique_target_path_list[i]))
    
    # Stage 2: extract the neighbourhood sequences
    seq_files, path_info_files = extract_target_neighbourhoods(
        params,
        params.input_gfa,
        target_seq_align_info,
        params.debug
    )

    # Stage 3: annotate the neighbourhoods
    all_seq_info_lists, _ = annotate_all_neighbourhoods(
        params, seq_files, path_info_files, unique_target_files, params.debug
    )

    # Stage 4: filter the annotations by coverage consistency
    trim_annotations_by_coverage(params, unique_target_files, all_seq_info_lists)

    LOG.info("Sarand ran successfully")
