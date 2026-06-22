"""
To be returned if a program version check returns an error
"""
PROGRAM_VERSION_NA = 'N/A'

"""
Python logging name for sarand
"""
SARAND_LOGGER_NAME = 'sarand'

"""
Global configuration for output file names / directories.

The output directory is laid out one sub-directory per pipeline stage:

    <output_dir>/
        target_hits/          (stage 1: target genes located in the graph)
        raw_neighborhoods/    (stage 2: extracted neighborhood sequences)
        final_neighborhoods/  (stages 3-4: annotated, coverage-filtered output)

Scratch/intermediate files live under the stage directory that produces them
and are removed at the end of the run unless --keep_intermediate_files is set.
"""

# Stage 1: target genes located in the assembly graph
TARGET_DIR_NAME = "target_hits"
TARGET_SEQ_DIR = "sequences"
TARGET_ALIGN_DIR = "alignments"
TARGET_ALIGN_FILE = "bandage.tsv"
TARGET_OVERLAP_FILE = "overlaps.txt"

# Stage 2: extracted neighborhood sequences
SEQ_DIR_NAME = "raw_neighborhoods"
NEIGHBORHOOD_SEQ_DIR = "neighborhood_sequences"
NEIGHBORHOOD_PATHS_DIR = "neighborhood_paths"
SEQ_NAME_PREFIX = "ng_sequences_"
# Intermediate files produced while enumerating/clustering neighborhood paths,
# grouped under <raw_neighborhoods>/intermediate_files/ (deleted by default).
INTERMEDIATE_DIR_NAME = "intermediate_files"
STREAM_PATHS_DIR = "stream_paths"             # per-direction (up/downstream) path clustering
MERGED_NEIGHBORHOOD_DIR = "merged_neighborhoods"  # full-neighborhood cd-hit inputs

# Stages 3-4: annotated, coverage-filtered neighborhoods
ANNOTATION_DIR = "final_neighborhoods"
