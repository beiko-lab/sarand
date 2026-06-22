"""Unit tests for sarand.external.cdhit command building."""
from pathlib import Path

from sarand.external.cdhit import CdhitParams


def test_cdhit_params_as_cmd():
    params = CdhitParams(
        input_file=Path("in.fasta"),
        output_file=Path("out.fasta"),
        identity=0.9,
        word_length=5,
        threads=4,
    )
    cmd = params.as_cmd()
    assert cmd[0] == "cd-hit"
    assert cmd[cmd.index("-c") + 1] == "0.9"
    assert cmd[cmd.index("-n") + 1] == "5"
    assert cmd[cmd.index("-T") + 1] == "4"
    # input/output passed as absolute paths
    assert cmd[cmd.index("-i") + 1] == str(Path("in.fasta").absolute())
    assert cmd[cmd.index("-o") + 1] == str(Path("out.fasta").absolute())
