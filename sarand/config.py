import os
from typing import Optional

"""
Configuration options for running dependencies through "conda run -n".

If any of these environment variables are supplied, then the respective
program will be run through "conda run -n name".
"""
CONDA_BAKTA_NAME: Optional[str] = os.environ.get('CONDA_BAKTA_NAME')
CONDA_RGI_NAME: Optional[str] = os.environ.get('CONDA_RGI_NAME')
CONDA_GRAPH_ALIGNER_NAME: Optional[str] = os.environ.get('CONDA_GRAPH_ALIGNER_NAME')
CONDA_BANDAGE_NAME: Optional[str] = os.environ.get('CONDA_BANDAGE_NAME')
CONDA_BLAST_NAME: Optional[str] = os.environ.get('CONDA_BLAST_NAME')
CONDA_EXE_NAME: Optional[str] = os.environ.get('CONDA_EXE_NAME', 'conda')
CONDA_BAKTA_DB: Optional[str] = os.environ.get('BAKTA_DB')

"""
To be returned if a program version check returns an error
"""
PROGRAM_VERSION_NA = 'N/A'

"""
Python logging name for sarand
"""
SARAND_LOGGER_NAME = 'sarand'

"""
Global configuration for output file names / directories
"""
AMR_DIR_NAME = "AMR_info"
AMR_SEQ_DIR = "sequences"
AMR_ALIGN_DIR = "alignments"
AMR_OVERLAP_FILE = "overlaps.txt"
SUBGRAPH_DIR_NAME = "subgraphs"
SEQ_DIR_NAME = "sequences_info"
SEQ_NAME_PREFIX = "ng_sequences_"
ANNOTATION_DIR = "annotations"
EVAL_DIR = "evaluation"
