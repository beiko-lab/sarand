#!/usr/bin/env python
import sys
import os
import argparse
import datetime
import logging
import shutil
import pkg_resources
from pathlib import Path

from sarand.__init__ import __version__
from sarand.full_pipeline import full_pipeline_main
from sarand.utils import (
    initialize_logger,
    check_dependencies,
    check_file,
    validate_range,
)


def main():
    """
    Main CLI entrypoint for sarand
    """
    run_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    parser = argparse.ArgumentParser(
        description="Identify and extract the "
        "local neighbourhood of target genes "
        " (such as AMR) from a GFA formatted "
        " assembly graph",
        prog="sarand",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "-i",
        "--input_gfa",
        required=True,
        help="Path to assembly graph (in GFA format) " "that you wish to analyse",
        type=check_file,
    )
    parser.add_argument(
        "-a",
        "--assembler",
        choices=["metaspades", "bcalm", "megahit"],
        required=True,
        help="Assembler used to generate input GFA (required to correctly parse "
        "coverage information)",
    )
    parser.add_argument(
        "-k",
        "--max_kmer_size",
        required=True,
        type=int,
        help="Maximum k-mer sized used by assembler to generate input GFA",
    )
    parser.add_argument(
        "--extraction_timeout",
        default=-1,
        type=int,
        help="Maximum time to extract neighbourhood, -1 indicates no limit",
    )
    parser.add_argument(
        "-j",
        "--num_cores",
        default=1,
        type=validate_range(int, 1, 100),
        help="Number of cores to use",
    )
    parser.add_argument(
        "-c",
        "--coverage_difference",
        default=30,
        type=validate_range(int, -1, 500),
        help="Maximum coverage difference to include "
        "when filtering graph neighbourhood. Use "
        "-1 to indicate no coverage threshold "
        "(although this will likely lead to false "
        "positive neighbourhoods).",
    )
    parser.add_argument(
        "-t",
        "--target_genes",
        default=Path(
            pkg_resources.resource_filename(__name__, "data/CARD_AMR_seq.fasta")
        ),
        type=check_file,
        help="Target genes to "
        "search for in the assembly graph (fasta formatted). "
        " Default is the pre-installed CARD database",
    )
    parser.add_argument(
        "-x",
        "--min_target_identity",
        default=95,
        type=validate_range(float, 0.1, 100),
        help="Minimum identity/coverage to identify presence "
        "of target gene in assembly graph",
    )
    parser.add_argument(
        "-l",
        "--neighbourhood_length",
        default=1000,
        type=validate_range(int, 0, 100000),
        help="Size of gene neighbourhood to extract from the assembly graph",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        help="Output folder for current run of sarand",
        default=Path(f"sarand_results_{run_time}"),
    )
    parser.add_argument(
        "-f",
        "--force",
        default=False,
        action="store_true",
        help="Force overwrite any previous files/output directories",
    )
    parser.add_argument(
        "--verbose",
        default=False,
        action='store_true',
        help='Provide verbose debugging output when logging',
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--no_rgi",
        default=False,
        action="store_false",
        help="Disable RGI based annotation of graph neighbourhoods",
    )
    group.add_argument(
        "--rgi_include_loose",
        default=False,
        action="store_true",
        help="Include loose criteria hits if using RGI to annotate"
        " graph neighbourhoods",
    )

    args = parser.parse_args()

    if os.path.exists(args.output_dir):
        if not args.force:
            parser.error(
                f"{args.output_dir} already exists, please use a different "
                "--output_dir or use --force to overwrite this directory"
            )
        else:
            print(f"Overwriting previously created {args.output_dir}", file=sys.stderr)
            shutil.rmtree(args.output_dir)
            os.makedirs(args.output_dir)
    else:
        os.makedirs(args.output_dir)

    initialize_logger(os.path.join(args.output_dir, f"run_{run_time}.log"), args.verbose)

    # check dependencies work
    # annoyingly Bandage --version requires X but Bandage --help does not
    dependencies = ["Bandage --help", "prokka --version", "blastn -version"]
    #cwd = os.getcwd()
    #PROKKA_COMMAND_PREFIX = 'docker run -v '+cwd+':/data staphb/prokka:latest '
    #dependencies = ["/media/Data/tools/Bandage_Ubuntu_dynamic_v0_8_1/Bandage --version",PROKKA_COMMAND_PREFIX+ "prokka --version", "blastn -version"]
    if not args.no_rgi:
        dependencies.append("rgi main --version")
    # check_dependencies(dependencies)


    # convert argparse to config dictionary
    args.run_time = run_time

    # logging file
    logging.info(f"Sarand initialised: output={args.output_dir}")

    # execute main workflow
    full_pipeline_main(args)


if __name__ == "__main__":
    main()
