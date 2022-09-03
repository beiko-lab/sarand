#!/bin/bash
set -euo pipefail
sarand -i test/spade_output/assembly_graph_with_scaffolds.gfa -o test/test_output -a metaspades -k 55

