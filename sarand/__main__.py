import argparse
import datetime
import shutil
import sys
from pathlib import Path

from sarand.__init__ import __version__
from sarand.pipeline import run_graph_pipeline
from sarand.util.logger import create_logger, get_logger
from sarand.util.cli import assert_dependencies_exist, check_file, validate_range
from sarand.databases import DATABASES, get_target_fasta, update_database

def main() -> None:
    """Main CLI entrypoint for sarand: parse arguments and dispatch a pipeline."""
    run_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    parser = argparse.ArgumentParser(
        description="Identify, extract, deduplicate, and coverage-filter the local "
                    "neighbourhoods of target genes (e.g., AMR genes) "
                    "from a GFA-formatted assembly graph",
        prog="sarand",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "-i",
        "--input_gfa",
        help="Path to assembly graph (in GFA format) that you wish to search",
        type=check_file,
    )
    parser.add_argument(
        "-a",
        "--assembler",
        choices=["metaspades", "bcalm", "megahit"],
        default=None,
        help="Assembler used to generate input GFA "
             "(required to correctly parse coverage information).",
    )
    parser.add_argument(
        "-k",
        "--max_kmer_size",
        type=int,
        help="Maximum k-mer sized used by assembler to generate input GFA "
             "(required to correctly calculate coverage).",
    )
    parser.add_argument(
        "--extraction_timeout",
        default=10000,
        type=int,
        help="Maximum number of minutes to spend traversing each target gene neighbourhood "
             " (high complexity subgraphs can be computationally demanding to traverse fully).",
    )
    parser.add_argument(
        "-j",
        "--num_cores",
        default=1,
        type=validate_range(int, 1, 100),
        help="Number of cores to use when running Sarand",
    )
    parser.add_argument(
        "-c",
        "--coverage_difference",
        default=30,
        type=validate_range(int, -1, 500),
        help="Maximum coverage difference within a path to retain when filtering graph neighbourhoods. "
             "Use -1 to indicate no coverage threshold (this will likely lead to chimeric false neighbourhoods)."
    )
    parser.add_argument(
        "-t",
        "--target_genes",
        default=None,
        type=check_file,
        help="Fasta-formatted nucleotide target gene sequences to search for in the assembly graph "
             "(Overrides --database or default CARD homolog sequences).",
    )
    parser.add_argument(
        "-d",
        "--database",
        choices=list(DATABASES),
        default="card",
        help="Reference target-gene database to search with when not supplying custom --target_genes."
    )
    parser.add_argument(
        "--update",
        action="store_true",
        default=False,
        help="Download/refresh the pre-installed CARD and NCBI AMR gene databases to their latest releases and exit."
    )
    parser.add_argument(
        "-x",
        "--min_target_identity",
        default=95,
        type=validate_range(float, 0.1, 100),
        help="Minimum identity for target gene hits in assembly graph",
    )

    parser.add_argument(
        "-c",
        "--min_target_coverage",
        default=95,
        type=validate_range(float, 0.1, 100),
        help="Minimum coverage for target gene hits in assembly graph",
    )
    parser.add_argument(
        "-l",
        "--neighbourhood_length",
        default=1000,
        type=validate_range(int, 0, 100000),
        help="Maximum gene neighbourhood length to extract surrounding each target gene hit in the assembly graph (bp).",
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
        "--deduplication_identity",
        default=0.9,
        type=validate_range(float, 0, 1),
        help="CD-HIT identity threshold for deduplicating extracted neighbourhoods",
    )


    # Parse arguments
    args = parser.parse_args()

    # --update just refreshes a reference database and exits; no graph needed.
    if args.update:
        create_logger(verbose=args.verbose)
        update_database(args.database)
        sys.exit(0)

    # Otherwise an assembler (and a graph / k-mer size) is required.
    if args.assembler is None:
        parser.error("one of the arguments -a/--assembler or --update is required")
    if args.max_kmer_size is None:
        parser.error("The --max_kmer_size (-k) argument is required.")

    # Override the keep intermediate files option if debug is set
    if args.debug:
        args.keep_intermediate_files = True

    # Setup the output logger path
    logger_output_path = Path(args.output_dir) / f"run_{run_time}.log"

    # Check if the output directory exists
    if Path(args.output_dir).exists():
        if not args.force:
            log = create_logger(verbose=args.verbose)
            log.error(
                f"{args.output_dir} already exists, please use a different "
                "--output_dir or use --force to overwrite this directory"
            )
            sys.exit(1)
        else:
            shutil.rmtree(args.output_dir)
            Path(args.output_dir).mkdir(parents=True)
            log = create_logger(output=logger_output_path, verbose=args.verbose)
            log.info(f"Overwriting previously created {args.output_dir}")

    else:
        Path(args.output_dir).mkdir(parents=True)
        create_logger(output=logger_output_path, verbose=args.verbose)

    # initialise the logger
    log = get_logger()

    # check dependencies work
    assert_dependencies_exist()

    # resolve the target genes: an explicit -t wins, otherwise use the selected
    # reference database (by default bundled CARD homolog nucleotide sequences).
    if args.target_genes is None:
        try:
            args.target_genes = get_target_fasta(args.database)
        except FileNotFoundError as e:
            log.error(str(e))
            sys.exit(1)
    log.info(f"Using target genes: {args.target_genes}")

    # insert run_time into config dictionary
    args.run_time = run_time

    # logging file
    log.info(f"Output directory: output={args.output_dir}")

    # execute the workflow for the chosen assembler
    run_graph_pipeline(args)


if __name__ == "__main__":
    main()
