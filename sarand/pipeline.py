"""Top-level orchestration of the main (assembly-graph) sarand pipeline.

Stages: find AMR genes in the graph -> extract their neighbourhoods -> ORF
annotate the neighbourhoods -> filter the annotations by coverage consistency.
"""
from __future__ import annotations

import argparse
import sys

from sarand.amr_finder import find_all_amr_in_graph
from sarand.annotation import annotate_all_amrs
from sarand.coverage import trim_annotations_by_coverage
from sarand.extract_neighborhood import sequence_neighborhood_main
from sarand.util.logger import LOG


def run_graph_pipeline(params: argparse.Namespace) -> None:
    """Run the full assembly-graph sarand workflow."""
    LOG.info("Starting analysis...")

    # Stage 1: find the target genes in the assembly graph
    LOG.info(f"Finding AMR genes in the assembly graph: {params.input_gfa}")
    unique_amr_files, unique_amr_path_list = find_all_amr_in_graph(
        params.input_gfa,
        params.output_dir,
        params.target_genes,
        params.min_target_identity,
        params.num_cores,
        params.keep_intermediate_files,
        params.debug,
    )

    if not unique_amr_files:
        LOG.error("No AMR genes were found in graph!")
        sys.exit(1)

    if not (unique_amr_path_list and len(unique_amr_path_list) == len(unique_amr_files)):
        LOG.error("AMR alignment info is not available")
        sys.exit(1)

    # pair each found AMR with its alignment info
    amr_seq_align_info = []
    for i, amr_file in enumerate(unique_amr_files):
        amr_seq_align_info.append((amr_file, unique_amr_path_list[i]))

    # Stage 2: extract the neighbourhood sequences
    seq_files, path_info_files = sequence_neighborhood_main(
        params,
        params.input_gfa,
        amr_seq_align_info,
        params.debug
    )

    # Stage 3: annotate the neighbourhoods
    all_seq_info_lists, _ = annotate_all_amrs(
        params, seq_files, path_info_files, unique_amr_files, params.debug
    )

    # Stage 4: filter the annotations by coverage consistency
    trim_annotations_by_coverage(params, unique_amr_files, all_seq_info_lists)

    LOG.info("Sarand ran successfully")
