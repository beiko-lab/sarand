#!/bin/bash
#
# Whole-pipeline regression test: run sarand on a bundled GFA and compare the
# key outputs of each stage against test/expected_output/ (a frozen run of the
# current code). Comparisons hash the sorted concatenation of each file group,
# so they are insensitive to file ordering and to the timestamps embedded in
# output filenames, but sensitive to any change in content.
#
# Regenerate the baseline deliberately (after an intended output change) with:
#   rm -rf test/expected_output && cp -r <a clean run> test/expected_output
#   find test/expected_output -name '*.clstr' -delete   # drop cd-hit intermediates
#
# Note: pipefail is intentionally NOT set so that a missing baseline file group
# hashes as empty and is reported as a failed comparison rather than aborting.
set -eu

GFA=test/spade_output/assembly_graph_with_scaffolds.gfa
OUT=test/actual_output
EXP=test/expected_output
NL=1000   # default --neighbourhood_length

# The cd-hit dedup step writes a temp_file.fasta into the working directory;
# make sure these strays are removed however the script exits.
trap 'rm -f temp_file.fasta temp_file.fasta.clstr' EXIT

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
compare "identified AMR genes"          "AMR_info/sequences/*.fasta"
compare "AMR overlap groups"            "AMR_info/overlaps.txt"
compare "extracted neighbourhoods"      "sequences_info/sequences_info_${NL}/sequences/ng_sequences_*.txt"
compare "neighbourhood path/coverage"   "sequences_info/sequences_info_${NL}/paths_info/ng_sequences_*.csv"
compare "ORF annotations"               "annotations/annotations_${NL}/*/annotation_detail_*.csv"
compare "unique ORF annotations"        "annotations/annotations_${NL}/*/trimmed_annotation_info_*.csv"
compare "coverage-filtered annotations" "annotations/annotations_${NL}/*/coverage_annotation_30_*.csv"

echo ""
if [[ "$fail" -ne 0 ]]; then
    echo "Functional test FAILED (output kept in $OUT for inspection)"
    exit 1
fi

rm -rf "$OUT"
echo "Functional test passed"
exit 0
