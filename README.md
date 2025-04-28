# Sarand

![sarand](sarand/docs/sarand.png)

Sarand is a tool to identify genes within an assembly graph and extract the local graph neighbourhood.
It has primarily been developed for the analysis of Antimicrobial Resistance (AMR) genes within metagenomic assembly graphs.
[CARD](card.mcmaster.ca) database is the default set of genes used for which neighborhoods are found but Sarand can support any user-supplied nucleotide fasta file of target genes.
<!--- Currently this is fixed to using the [CARD](card.mcmaster.ca) database but will be expanded in the near future to support any user-supplied nucleotide fasta file of target genes.-->


![sarand overview](sarand/docs/sarand_summary.png)

## 1. Installation

Sarand can be run using a conda environment or in a container (Docker or Singularity) and requires 4 key dependencies:

- [Bakta](https://github.com/oschwengers/bakta)
- [RGI](https://github.com/arpcard/rgi)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Bandage](https://github.com/rrwick/Bandage) or [GraphAligner](https://github.com/maickrau/GraphAligner)

### 1a. Docker

This is the easiest way to run Sarand, note that the `-v` argument maps a host directory to the Docker container.
You need to replace `/host/path` and `/container/path` in the command below with the path to the directory containing your input GFA.
Note that this will also be the location that the output is written to.

The most simple way to approach this is by mapping `/host/path` and `/container/path` to the same directory to keep paths consistent.

```shell
docker run -v /host/path:/container/path -it beiko-lab/sarand:1.1.1 -i /container/path/input.gfa -o /container/path/output
```

### 1b. Singularity

As singularity will automatically map paths, you simply need to run it in the format of:

```shell
singularity run docker://beiko-lab/sarand:1.1.1 -i input.gfa -o output
```


### 1c. Conda

As there are dependency conflicts between the tools used by sarand, you will need to create multiple conda environments.

**Creating environments:**

```shell
# 1. Create the sarand environment
conda create -n sarand-1.1.1 -c conda-forge -c bioconda -y blast=2.14.0 dna_features_viewer=3.1.2 numpy matplotlib-base gfapy=1.2.3 cd-hit=4.6.8 networkx gzip pandas python pillow biopython

# 2. Create the bakta environment
conda create -n bakta-1.8.1 -c conda-forge -c bioconda -y bakta=1.8.1

# 3.a. Create the Bandage environment
conda create -n bandage-0.8.1 -c conda-forge -c bioconda -c defaults -y bandage=0.8.1

# 3.b. Create the GraphAligner environment
conda create -n graphaligner-1.0.17b -c conda-forge -c bioconda -y graphaligner=1.0.17b

# 4. Create the RGI environment
conda create -n rgi-5.2.0 -c conda-forge -c bioconda -c defaults -y rgi=5.2.0
```
Please note that step 3.b is not required if you run the default version of Sarand. Sarand, by default, utilizes Bandage for sequence alignment in the assembly graphs. However, If you prefer to use GraphAligner, please make sure to run command 3.b and install it.

**Downloading and updating the Bakta database:**

```shell
cd /tmp
wget https://zenodo.org/record/7669534/files/db-light.tar.gz
tar -xzvf db-light.tar.gz
rm db-light.tar.gz

# Note: Here you will need to specify a path to keep the Bakta database
# This example uses /db/bakta but you can use any path you like
mkdir -p /db/bakta
mv db-light /db/bakta
conda run -n bakta-1.8.1 amrfinder_update --force_update --database /db/bakta/db-light/amrfinderplus-db
```

**Configuring conda environments:**

Here you will specify environment variables that are specific to the `sarand-1.1.1` environment,
these will be automatically used when the environment is active.

```shell
conda activate sarand-1.1.1
conda env config vars set CONDA_BAKTA_NAME=bakta-1.8.1
conda env config vars set CONDA_BANDAGE_NAME=bandage-0.8.1
conda env config vars set CONDA_RGI_NAME=rgi-5.2.0
conda env config vars set BAKTA_DB=/db/bakta/db-light
# Note1: Only run the following command if you have created graphaligner-1.0.17b conda environemnt in step 3.b above.
conda env config vars set CONDA_GRAPH_ALIGNER_NAME=graphaligner-1.0.17b

# Note2: Here you can specify an alternate exe (e.g. micromamba, mamba).
conda env config vars set CONDA_EXE_NAME=conda
```

**Installing sarand:**

```shell
conda activate sarand-1.1.1
# python -m pip install sarand==1.1.1
pip install sarand==1.1.1
```

## 2. Testing

You can test your install has worked by running the test script via `bash test/test.sh`
This will execute sarand on a test dataset (using the following command) and check all the expected outputs are created correctly.

    `sarand -i test/spade_output/assembly_graph_with_scaffolds.gfa -o test/test_output -a metaspades -k 55`



## 3. Usage

All of sarand's parameters can be set using the command line flags.
The only required input file is an assembly graph in `.gfa` format.

This can be generated using metagenomic (or genomic) de-novo assembly tools
such as [metaSPAdes](https://github.com/ablab/spades) or [megahit](https://github.com/voutcn/megahit).
If your chosen assembly tool generates a `fastg` formatted graph utilities such as `fastg2gfa` can be used to convert them.

```
usage: sarand [-h] [-v] -i INPUT_GFA -a ASSEMBLER
              -k MAX_KMER_SIZE [-j NUM_CORES] [-c COVERAGE_DIFFERENCE]
              [-t TARGET_GENES] [-x MIN_TARGET_IDENTITY]
              [-l NEIGHBOURHOOD_LENGTH] [-o OUTPUT_DIR] [-f]
              [--verbose] [--no_rgi | --rgi_include_loose] [--use_ga]
              [--ga] [--keep_intermediate_files] [--debug]

Identify and extract the local neighbourhood of target genes (such as AMR)
from a GFA formatted assembly graph

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i INPUT_GFA, --input_gfa INPUT_GFA
                      Path to assembly graph (in GFA format) that you wish
                      to analyse
  -a {metaspades,bcalm,megahit,metacherchant,contig},
  --assembler {metaspades,bcalm,megahit,metacherchant,contig}
                      Assembler used to generate input GFA (required to
                      correctly parse coverage information)
  -k MAX_KMER_SIZE, --max_kmer_size MAX_KMER_SIZE
                      Maximum k-mer sized used by assembler to generate
                      input GFA
  --extraction_timeout EXTRACTION_TIMEOUT
                      Maximum time to extract neighbourhood per gene in
                      minutes, -1 indicates no limit
  -j NUM_CORES, --num_cores NUM_CORES
                      Number of cores to use
  -c COVERAGE_DIFFERENCE, --coverage_difference COVERAGE_DIFFERENCE
                      Maximum coverage difference to include when filtering
                      graph neighbourhood. Use -1 to indicate no coverage
                      threshold (although this will likely lead to false
                      positive neighbourhoods).
  -t TARGET_GENES, --target_genes TARGET_GENES
                      Target genes to search for in the assembly graph
                      (fasta formatted). Default is the pre-installed CARD database
  -x MIN_TARGET_IDENTITY, --min_target_identity MIN_TARGET_IDENTITY
                      Minimum identity/coverage to identify presence of
                      target gene in assembly graph
  -l NEIGHBOURHOOD_LENGTH, --neighbourhood_length NEIGHBOURHOOD_LENGTH
                      Size of gene neighbourhood to extract from the
                      assembly graph
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                      Output folder for current run of sarand
  -f, --force         Force overwrite any previous files/output     
                      directories
  --verbose           Provide verbose debugging output when logging,
                      and keep intermediate files
  --no_rgi            Disable RGI based annotation of graph neighbourhoods
  --rgi_include_loose Include loose criteria hits if using RGI to annotate
                      graph neighbourhoods
  --use_ga            Enable GraphAligner (instead of Bandage) for  
                      sequence alignment in the graph
  --ga [GA ...]       Additional arguments to supply to graph aligner in the
                      form of --ga key value, e.g. --ga E-cutoff 0.1;
                      it should be used only if use_ga is set to True
  --keep_intermediate_files
                      Do not delete intermediate files.
  --debug               Creates additional files for debugging purposes.
  -seq SEQ_NUMBER, --max_number_seq_for_cdhit SEQ_NUMBER    
  		      Max Number of sequence for cd-hit
  -sim [0 1],  -similarity [0 1]
                     similarity threshold for cdhit (a number between 0 and 1)
  --meta_main_dir METACHERCHANT_MAIN_DIR
  		     The main directory for metacherchant containing
  		     AMR_seqs_full.fasta, all AMR sequences and the
  		     extracted local graphs by metacherchant.
```

**Running for Metacherchant:**

To extract neighborhoods from Metacherchant, you first need to run Metacherchant separately on your set of target genes. For each gene, Metacherchant will generate a local neighborhood graph. However, it does not provide the actual neighborhood sequences. To extract these sequences from the generated local graphs, Sarand must be run on them.

***Required Input Files and Directories***

Ensure that the following items are placed inside a directory, which will be passed to Sarand as `meta_main_dir`:
* `AMR_seqs_full.fasta`: a FASTA file containing all target genes' names and sequences (each in a separate line)
* `AMR_info/sequences`: a directory containing separate files each each named <geneName>.fasta, where each file contains the name and sequence of a target gene.
* `output` a directory generated by Metacherchant containing the local graphs produced for the genes.

***Running Sarand***
Once the required files and directories are set up, execute Sarand for Metacherchant using the following command:
```shell
sarand -o <output_dir> -a metacherchant --meta_main_dir <meta_main_dir>
```

**Running for Contigs:**

To extract neighborhoods from contigs, execute Sarand using the following command:
```shell
sarand -i <contigs_file> -o <output_dir> -a contig
```

### 3a. Output
All results will be available in specified output directory (default is `sarand_results_` followed by a timestamp).
Here is the list of important directories and files that can be seen there and a short description of their content:
* `AMR_info`: this directory contains the list of identified AMR sequences.
    * `AMR_info/sequences/`:The sequence of identified AMRs, from graph, is stored here, with a name similar to their original name (file name is generated by calling `sarand/utils.py::restricted_amr_name_from_modified_name(amr_name_from_title(amr_original_name)))`
    * `AMR_info/alignments/`: The alignment details for all AMR sequences are stored here.

* `sequences_info/sequences_info_{neighbourhood_length}/`: This directory stores the information of extracted neighborhood sequences from the assembly graph.
    * `sequences_info/sequences_info_{params.neighbourhood_length}/sequences/`: the extracted sequences in the neighborhood of each AMR are stored in a file like `ng_sequences_{AMR_NAME}_{params.neighbourhood_length}_{DATE}.txt`.
For each extracted sequence, the first line denotes the corresponding path, where the nodes representing the AMR sequence are placed in '[]'. The next line denotes the extracted sequence where the AMR sequence is in lower case letters and the neighborhood is in upper case letters.
    * `sequences_info/sequences_info_{params.neighbourhood_length}/paths_info/`: The information of nodes representing the AMR neighborhood including their name, the part of the sequence represented by each node (start position and end position) as well as their coverage is stored in a file like `ng_sequences_{AMR_NAME}_{params.neighbourhood_length}_{DATE}.csv`

* `annotations/annotations_{params.neighbourhood_length}`: The annotation details are stored in this directory.
    * `annotations/annotations_{params.neighbourhood_length}/annotation_{AMR_NAME}_{params.neighbourhood_length}`: this directory contains all annotation details for a given AMR.
    * `gene_comparison_<AMR_NAME>.png`: An image visualizing annotations
    * `annotation_detail_{AMR_NAME}.csv`: the list of annotations of all extracted sequences for an AMR gene
    * `trimmed_annotation_info_{AMR_NAME}.csv`: the list of unique annotations of all extracted sequences for an AMR gene
    * `coverage_annotation_{COVERAGE_DIFFERENCE}_{AMR_NAME}.csv`: the list of the annotations in which the gene coverage difference from the AMR gene coverage is less than GENE_COVERAGE_DIFFERENCE value.
    * `bakta_dir_extracted{NUM}_{DATE}`: it contains the output of prokka for annotation of a sequence extracted from the neighborhood of the target AMR gene in the assembly graph.
    * `rgi_dir`: contains RGI annotation details for all extracted neighborhood sequences of the target AMR gene.
