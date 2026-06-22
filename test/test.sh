#!/bin/bash
#
# Integrated Test driver for sarand. Runs in two stages:
#   1. Unit tests (pytest test/unit) covering the pure functions of each module.
#   2. A whole-pipeline regression test: run sarand on a bundled GFA and compare
#      the key outputs of each stage against test/expected_output/ (a frozen run
#      of the current code). Comparisons hash the sorted concatenation of each
#      file group, so they are insensitive to file ordering and to the
#      timestamps embedded in output filenames, but sensitive to any change in
#      content.
#
# The unit tests need pytest (pip install -e '.[test]'); the functional test
# additionally needs the sarand binary and its external tools on the PATH.
#
# Regenerate the baseline deliberately (after an intended output change) with:
#   rm -rf test/expected_output && cp -r <a clean run> test/expected_output
#   find test/expected_output -name '*.clstr' -delete   # drop cd-hit intermediates
#   find test/expected_output -name '*.log'   -delete   # drop per-run logs
#   rm -rf test/expected_output/target_hits/alignments  # drop bandage intermediate
#
# Note: not using pipefail as a missing baseline file group will cause an abortion
# this way it will just hash empty
set -eu

GFA=test/spade_output/assembly_graph_with_scaffolds.gfa
OUT=test/actual_output
EXP=test/expected_output
# default --neighborhood_length
NL=1000   

# The cd-hit dedup step writes a temp_file.fasta into the working directory;
# make sure these strays are removed however the script exits.
trap 'rm -f temp_file.fasta temp_file.fasta.clstr' EXIT

echo ""
echo "Running unit tests ..."
python -m pytest test/unit -q

echo ""
echo "Running functional/integration test ..."
rm -rf "$OUT"
sarand -i "$GFA" -o "$OUT" -a metaspades -k 55

fail=0
compare() {
    # $1 = human-readable label, $2 = file glob relative to the output root
    local label="$1" glob="$2" exp act
    exp=$(cat $EXP/$glob 2>/dev/null | sort | md5sum)
    act=$(cat $OUT/$glob 2>/dev/null | sort | md5sum)
    if [[ "$exp" == "$act" ]]; then
        echo "  ok   - $label"
    else
        echo "  FAIL - $label  (glob: $glob)"
        fail=1
    fi
}

echo ""
echo "Comparing outputs against $EXP ..."
compare "identified target genes"       "target_hits/sequences/*.fasta"
compare "target overlap groups"         "target_hits/overlaps.txt"
compare "extracted neighborhoods"       "raw_neighborhoods/neighborhood_sequences/ng_sequences_*.txt"
compare "neighborhood path/coverage"    "raw_neighborhoods/neighborhood_paths/ng_sequences_*.csv"
compare "ORF annotations"               "final_neighborhoods/annotation_*_${NL}/annotation_detail_*.csv"
compare "called ORFs (gff)"             "final_neighborhoods/annotation_*_${NL}/orfs_*.gff"
compare "coverage-filtered annotations" "final_neighborhoods/annotation_*_${NL}/coverage_annotation_30_*.csv"
compare "combined final neighborhoods"  "final_neighborhoods/final_neighborhoods.csv"

echo ""
if [[ "$fail" -ne 0 ]]; then
    echo "Functional test FAILED (output kept in $OUT for inspection)"
    exit 1
fi

rm -rf "$OUT"
echo "Functional test passed"
exit 0
