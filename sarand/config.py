import os
import socket

"""
This will check if sarand is running in a docker container, this is set in 
the dockerfile. If it is running in a docker container, then subprocess
calls will use the conda environment specific to each dependency as specified.
"""

# True if the script is being run in the sarand docker container
IS_DOCKER_CONTAINER: bool = os.environ.get('IS_DOCKER_CONTAINER', 'False') == '1'

DOCKER_BAKTA_ENV = 'bakta'
DOCKER_BAKTA_ENV = 'bakta-1.8.1'

DOCKER_RGI_ENV = 'rgi'
DOCKER_RGI_ENV = 'rgi-5.2.0'
# DOCKER_RGI_ENV = 'rgi-6.0.2'

DOCKER_GRAPH_ALIGNER_ENV = 'graphaligner'
DOCKER_GRAPH_ALIGNER_ENV = 'graphaligner-1.0.17b'


DOCKER_CONDA_EXE = 'micromamba'
DOCKER_CONDA_EXE = 'conda'

# TODO: Expects BAKTA_DB to be provided as an environment variable


TMP_PROKKA_ENV = 'prokka_1.14.6'
TMP_BANDAGE_ENV = 'bandage-0.8.1'

DEBUG = False

IS_LAPTOP = socket.gethostname() == 'Aarons-MacBook-Pro.local'

if IS_LAPTOP:
    IS_DOCKER_CONTAINER = True
    DOCKER_BAKTA_ENV = 'bakta-1.8.1'
    DOCKER_RGI_ENV = 'rgi-5.2.0'
    DOCKER_GRAPH_ALIGNER_ENV = 'graphaligner-1.0.17b'
    DOCKER_CONDA_EXE = 'conda'
    TMP_PROKKA_ENV = 'prokka-1.14.5'
    TMP_BANDAGE_ENV = 'bandage-0.8.1'


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