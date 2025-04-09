import argparse
import datetime
import os
import shutil
import sys
from pathlib import Path

from sarand.__init__ import __version__
from sarand.full_pipeline import full_pipeline_main
from sarand.util.logger import create_logger, get_logger
from sarand.util.pkg import get_pkg_card_fasta_path
from sarand.utils import assert_dependencies_exist, check_file, validate_range
from sarand.test_metacherchant import test_metacherchant_main
from sarand.amr_neighborhood_in_contigs import find_contig_amrs_main

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
        #required=True,
        help="Path to assembly graph (in GFA format) " "that you wish to analyse",
        type=check_file,
    )
    parser.add_argument(
        "-a",
        "--assembler",
        choices=["metaspades", "bcalm", "megahit", "metacherchant", "contig"],
        required=True,
        help="Assembler used to generate input GFA (required to correctly parse "
             "coverage information)",
    )
    parser.add_argument(
        "-k",
        "--max_kmer_size",
        #required=True,
        type=int,
        help="Maximum k-mer sized used by assembler to generate input GFA",
    )
    parser.add_argument(
        "--extraction_timeout",
        default=10000,
        type=int,
        help="Maximum time to extract neighbourhood per AMR gene in minutes",
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
        default=get_pkg_card_fasta_path(),
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
        help='Provide verbose debugging output when logging, and keep intermediate files',
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
    parser.add_argument(
        '--use_ga',
        default = False,
        action='store_true',
        help="Enable GraphAligner (instead of Bandage) for sequence alignment in the graph",
    )
    # GraphAligner options
    parser.add_argument(
        '--ga',
        default=None,
        action='append',
        nargs='*',
        help='Additional arguments to supply to graph aligner in the form of --ga key value,'
        'e.g. --ga E-cutoff 0.1; it should be used only if use_ga is set to True'
    )
    parser.add_argument(
        "--keep_intermediate_files",
        default=False,
        action="store_true",
        help="Do not delete intermediate files.",
    )
    parser.add_argument(
        '--debug',
        default=False,
        action='store_true',
        help='Creates additional files for debugging purposes.',
    )
    parser.add_argument(
        "-n",
        "--max_down_up_paths",
        default=-1,
        type=validate_range(int, 0, 100000),
        help="Max Number of downstream or upstream paths",
    )

    parser.add_argument(
        "-seq",
        "--max_number_seq_for_cdhit",
        default=100000,
        type=validate_range(int, 0, 100000),
        help="Max Number of sequence for cd-hit",
    )

    parser.add_argument(
        "-sim",
        "--similarity",
        default=0.9,
        type=validate_range(float, 0, 1),
        help="similarity threshold for cdhit (a number between 0 and 1)",
    )
    parser.add_argument(
        '--meta_main_dir',
        default = None,
        help = 'The main directory for metacherchant containing AMR_seqs_full.fasta,'
            ' all AMR sequences and the extracted local graphs by metacherchant'
    )

    # Parse arguments
    args = parser.parse_args()

    # Enforce conditional requirement for input_gfa and max_kmer_size
    if args.assembler != "metacherchant" and args.input_gfa is None:
        parser.error("The --input_gfa (-i) argument is required.")
    if args.assembler != "metacherchant" and args.assembler != "contig" and args.max_kmer_size is None:
        parser.error("The --max_kmer_size (-k) argument is required.")
        
    # Override the keep intermediate files option if debug is set
    if args.debug:
        args.keep_intermediate_files = True

    # Setup the output logger path
    logger_output_path = os.path.join(args.output_dir, f"run_{run_time}.log")

    # Check if the output directory exists
    if os.path.exists(args.output_dir):
        if not args.force:
            log = create_logger(verbose=args.verbose)
            log.error(
                f"{args.output_dir} already exists, please use a different "
                "--output_dir or use --force to overwrite this directory"
            )
            sys.exit(1)
        else:
            shutil.rmtree(args.output_dir)
            os.makedirs(args.output_dir)
            log = create_logger(output=logger_output_path, verbose=args.verbose)
            log.info(f"Overwriting previously created {args.output_dir}")

    else:
        os.makedirs(args.output_dir)
        create_logger(output=logger_output_path, verbose=args.verbose)

    # Get the logger
    log = get_logger()

    # check dependencies work
    assert_dependencies_exist(graph_aligner = args.use_ga,
            bandage = not args.use_ga, rgi=not args.no_rgi)

    # convert argparse to config dictionary
    args.run_time = run_time

    # logging file
    log.info(f"Sarand initialized: output={args.output_dir}")

    # execute main workflow
    #full_pipeline_main(args)
    if args.assembler == "metacherchant":
        test_metacherchant_main(args)
    elif args.assembler == "contig":
        find_contig_amrs_main(args)
    else:
        # execute main workflow
        full_pipeline_main(args)


if __name__ == "__main__":
    main()
