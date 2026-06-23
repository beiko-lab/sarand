# Sarand


[![DOI](https://zenodo.org/badge/482245454.svg)](https://doi.org/10.5281/zenodo.20394760)


![sarand](sarand/docs/sarand.png)

Sarand is a tool to identify genes within an assembly graph and extract the local graph neighborhood.
It has primarily been developed for the analysis of Antimicrobial Resistance (AMR) genes within metagenomic assembly graphs.
[CARD](card.mcmaster.ca) or [NCBI Reference Gene Catalog](https://www.ncbi.nlm.nih.gov/pathogens/refgene/#) databases are built-in for searching in a graph but Sarand can support any user-supplied nucleotide fasta file of target genes.

![sarand overview](sarand/docs/sarand_summary.png)

## 1. Installation

Sarand can be run using a conda environment or in a container (Docker or Singularity) and requires these key external dependencies:

- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (not called directly but used by Bandage's graph alignment)
- [Bandage](https://github.com/rrwick/Bandage)
- [minimap2](https://github.com/lh3/minimap2) 
- [cd-hit](https://github.com/weizhongli/cdhit)


### 1a. Conda

All of sarand's dependencies are installed into a single conda environment;
Bandage and BLAST+ are simply expected on the `PATH`.

**Install dependencies:**

```shell
conda env create -f conda_env.yml
```

Which can alternative be run explicitly without using the enviornment yaml like this:

```shell
conda create -n sarand -c conda-forge -c bioconda -y blast=2.17.0 bandage=0.9.0 minimap2=2.31 gfapy=1.2.3 cd-hit=4.8.1 networkx biopython pyrodigal
```

**Install sarand:**

```shell
git clone https://github.com/beiko-lab/sarand.git
cd sarand
conda activate sarand
python -m pip install .
```

### 1b. Apptainer/Singularity

This is the easiest way to run Sarand. As apptainer/singularity will automatically map paths, you simply need to run it in the format of:

```shell
singularity run docker://ghcr.io/beiko-lab/sarand:latest -i input.gfa -o output -a metaspades -k 55
```

### 1c. Docker

A docker container is also provided but note you'll need to provide an argument to the `-v` option in docker as this maps a host directory to the Docker container.
Basically, just replace `/host/path` and `/container/path` in the command below with the path to the directory containing your input GFA.
Note that this will also be the location that the output is written to.

```shell
docker run -v /host/path:/container/path -it ghcr.io/beiko-lab/sarand:latest -i /container/path/input.gfa -o /container/path/output -a metaspades -k 55
```

## 2. Testing

You can test your install has worked by running the test script via `bash test/test.sh`
This will execute the unit tests and run sarand on a test dataset (using the following command) and check all the expected outputs are created correctly.

    sarand -i test/minimal_graph.gfa -o test/actual_output -a metaspades -k 55


## 3. Usage

All of sarand's parameters can be set using the command line flags.
The only required input file is an assembly graph in `.gfa` format.

This can be generated using metagenomic (or genomic) de-novo assembly tools
such as [metaSPAdes](https://github.com/ablab/spades) or [megahit](https://github.com/voutcn/megahit).
If your chosen assembly tool generates a `fastg` formatted graph utilities such as `fastg2gfa` can be used to convert them.

```
> sarand --help

usage: sarand [-h] [-v] [-i INPUT_GFA] [-a {metaspades,bcalm,megahit}] [-k MAX_KMER_SIZE]
              [--extraction_timeout EXTRACTION_TIMEOUT] [-j NUM_CORES]
              [-c COVERAGE_DIFFERENCE] [-t TARGET_GENES] [-d {card,ncbi}] [--update]
              [-x MIN_TARGET_IDENTITY] [-z MIN_TARGET_COVERAGE] [-l NEIGHBORHOOD_LENGTH]
              [-o OUTPUT_DIR] [-f] [--verbose] [--keep_intermediate_files] [--debug]
              [--deduplication_identity DEDUPLICATION_IDENTITY]

Identify, extract, deduplicate, and coverage-filter the local neighborhoods of target
genes (e.g., AMR genes) from a GFA-formatted assembly graph

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i, --input_gfa INPUT_GFA
                        Path to assembly graph (in GFA format) that you wish to search
  -a, --assembler {metaspades,bcalm,megahit}
                        Assembler used to generate input GFA (required to correctly parse
                        coverage information).
  -k, --max_kmer_size MAX_KMER_SIZE
                        Maximum k-mer sized used by assembler to generate input GFA
                        (required to correctly calculate coverage).
  --extraction_timeout EXTRACTION_TIMEOUT
                        Maximum number of minutes to spend traversing each target gene
                        neighborhood (high complexity subgraphs can be computationally
                        demanding to traverse fully).
  -j, --num_cores NUM_CORES
                        Number of cores to use when running Sarand
  -c, --coverage_difference COVERAGE_DIFFERENCE
                        Maximum coverage difference within a path to retain when filtering
                        graph neighborhoods. Use -1 to indicate no coverage threshold
                        (this will likely lead to chimeric false neighborhoods).
  -t, --target_genes TARGET_GENES
                        Fasta-formatted nucleotide target gene sequences to search for in
                        the assembly graph (Overrides --database or default CARD homolog
                        sequences).
  -d, --database {card,ncbi}
                        Reference target-gene database to search with when not supplying
                        custom --target_genes.
  --update              Download/refresh the pre-installed CARD and NCBI AMR gene
                        databases to their latest releases and exit.
  -x, --min_target_identity MIN_TARGET_IDENTITY
                        Minimum identity for target gene hits in assembly graph
  -z, --min_target_coverage MIN_TARGET_COVERAGE
                        Minimum coverage for target gene hits in assembly graph
  -l, --neighborhood_length NEIGHBORHOOD_LENGTH
                        Maximum gene neighborhood length radius (i.e., the number of bases
                        upstream and/or downstream) to extract surrounding each target
                        gene hit in the assembly graph (bp).
  -o, --output_dir OUTPUT_DIR
                        Output folder for current run of sarand
  -f, --force           Force overwrite any previous files/output directories
  --verbose             Provide verbose debugging output when logging, and keep
                        intermediate files
  --keep_intermediate_files
                        Do not delete intermediate files.
  --debug               Enable debug-level logging and create additional files for
                        debugging purposes.
  --deduplication_identity DEDUPLICATION_IDENTITY
                        CD-HIT identity threshold for deduplicating extracted
                        neighborhoods
```


**Reference databases (`--database` / `--update`):**

By default Sarand searches for the CARD antimicrobial-resistance genes ("protein homolog models").  Two reference databases of nucleotide target genes are
supported and can be kept up to date (other datasets can be used by supplying a fasta with `--target_genes`):

* `card` (default) : the [CARD](https://card.mcmaster.ca) protein-homolog model
  nucleotide sequences.
* `ncbi` : the [NCBI AMRFinderPlus](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/)
  reference gene CDS (`AMR_CDS.fa`).

`--update` downloads the latest release of the bundled databases (only if newer
than the local copy) and exits; the databases are cached in a user-writable
directory (`$SARAND_DB_DIR`, else `$XDG_DATA_HOME/sarand`, else
`~/.local/share/sarand`) and used automatically on subsequent runs.

```shell
# refresh the bundled CARD database to the latest release
sarand --update
sarand -i input.gfa -o output -a metaspades -k 55 --database ncbi
```

### 3a. Output
All results are written to the specified output directory (default is `sarand_results_`
followed by a timestamp). `{TARGET}` below is the (filesystem-safe) target gene name and
`{L}` is the `--neighborhood_length`.

* `target_hits/`: the target genes located in the assembly graph.
    * `target_hits/sequences/`: the target gene sequences recovered from the graph, one FASTA per target.
    * `target_hits/overlaps.txt`: groups of target genes whose graph paths overlap.
    * `target_hits/alignments/`: the Bandage alignment details (only with `--debug` or `--keep_intermediate_files`).

* `raw_neighborhoods/`: the neighborhood sequences extracted from the assembly graph.
    * `raw_neighborhoods/neighborhood_sequences/`: the extracted sequences for each target in a FASTA file like `ng_sequences_{TARGET}_{L}_{DATE}.fasta`. For each sequence, the first line is the node path (the nodes representing the target gene are wrapped in `()`) and the second line is the sequence, with the target gene in lower case and the flanking neighborhood in upper case.
    * `raw_neighborhoods/neighborhood_paths/`: the per-node path/coverage info backing each extracted sequence (node name, the start/end offset it covers, its coverage, and the whole-path coverage) in a file like `ng_sequences_{TARGET}_{L}_{DATE}.csv`.

* `final_neighborhoods/`: the annotated, coverage-filtered neighborhoods. The neighborhoods are ORF-annotated with pyrodigal, which calls open reading frames but does not assign gene names or functional products; the called ORFs are therefore reported directly as the `orfs_{TARGET}.{ffn,faa,gff}` files below rather than as a named-gene table.
    * `final_neighborhoods.fasta`: a single combined FASTA of every final neighborhood sequence across all targets (header `>{TARGET}_{seq_name}`).
    * `final_neighborhoods.csv`: a single combined summary, one row per final neighborhood, with columns `target_name, seq_name, target_gene, gene_path, target_coverage, coverages` (where `target_gene` is the target ORF's `start-end` coordinates and `gene_path` lists every ORF's `start-end`).
    * `final_neighborhoods/annotation_{TARGET}_{L}/`: the per-target details:
        * `orfs_{TARGET}.ffn` / `.faa` / `.gff`: the pyrodigal ORFs of the target's final neighborhoods as nucleotide FASTA, protein FASTA and GFF3.
        * `coverage_annotation_{COVERAGE_DIFFERENCE}_{TARGET}.csv`: the per-ORF coverage table of the neighborhoods kept after coverage-consistency filtering (only when `--coverage_difference > 0`); near-identical neighborhoods are collapsed.
    * `not_found_annotation_targets_in_graph.txt`: targets/sequences for which no annotation could be produced.

<!--
### 3b. Visualising annotations

Sarand no longer renders annotation comparison images as part of the main run.
A standalone helper script is provided to generate them on demand from any of
the annotation CSVs above (for example
`coverage_annotation_{COVERAGE_DIFFERENCE}_{AMR_NAME}.csv`):

```shell
python scripts/visualize_annotation.py --csvfile <annotation.csv> --output gene_comparison.png --title "<AMR_NAME>"
```

The script is not installed with the `sarand` package; run it directly from a
checkout. It only needs `dna_features_viewer`, `matplotlib`, `numpy` and
`pillow`, which are already dependencies of sarand.
-->
