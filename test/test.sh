#!/bin/bash
set -euo pipefail
sarand -i test/spade_output/assembly_graph_with_scaffolds.gfa -o test/actual_output -a metaspades -k 55 

# check alignments match
expected=$(cat test/expected_output/AMR_info/alignments/*.tsv | sort | md5sum)
actual=$(cat test/actual_output/AMR_info/alignments/*.tsv | sort | md5sum)
if [[ "$expected" != "$actual" ]]
then
    echo "Functional test 1: error in bandage alignments"
    exit 1
fi

expected=$(cat test/expected_output/sequences_info/sequences_info_1000/sequences/*.txt | sort | md5sum)
actual=$(cat test/actual_output/sequences_info/sequences_info_1000/sequences/*.txt | sort | md5sum)
if [[ "$expected" != "$actual" ]]
then
    echo "Functional test 2 failed: mismatch in extracted sequences"
    exit 1
fi

expected=$(cat test/expected_output/annotations/annotations_1000/*/*.csv | sort | md5sum)
actual=$(cat test/actual_output/annotations/annotations_1000/*/*.csv | sort | md5sum)
if [[ "$expected" != "$actual" ]]
then
    echo "Functional test 3 failed: mismatch in annotations"
    exit 1
fi

rm -rf test/actual_output

echo -e "\n\nFunctional tests passed: 3/3\n\n"
exit 0
