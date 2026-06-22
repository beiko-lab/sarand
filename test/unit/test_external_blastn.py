"""Unit tests for sarand.external.blastn parsing and command building."""
from pathlib import Path

import pytest

from sarand.external.blastn import Blastn, BlastnOutFmt, BlastnParams, BlastnResult


def test_blastn_params_as_cmd_minimal():
    params = BlastnParams(
        query=Path("q.fasta"),
        subject=Path("s.fasta"),
        task="blastn-short",
        outfmt=BlastnOutFmt.FMT_1,
    )
    cmd = params.as_cmd()
    assert cmd[0] == "blastn"
    assert "-task" in cmd and "blastn-short" in cmd
    assert "-outfmt" in cmd
    assert cmd[cmd.index("-outfmt") + 1] == "10 pident length qcovhsp"
    # optional filters absent
    assert "-evalue" not in cmd
    assert "-perc_identity" not in cmd


def test_blastn_params_as_cmd_with_optional_filters():
    params = BlastnParams(
        query=Path("q.fasta"),
        subject=Path("s.fasta"),
        task="blastn",
        outfmt=BlastnOutFmt.FMT_1,
        max_target_seqs=5,
        evalue=1e-10,
        perc_identity=95.0,
    )
    cmd = params.as_cmd()
    assert cmd[cmd.index("-evalue") + 1] == "1e-10"
    assert cmd[cmd.index("-perc_identity") + 1] == "95.0"
    assert cmd[cmd.index("-max_target_seqs") + 1] == "5"


def test_blastn_result_from_outfmt_fmt1():
    res = BlastnResult.from_outfmt("99.5,150,98", BlastnOutFmt.FMT_1)
    assert res.pident == 99.5
    assert res.length == 150
    assert res.qcovhsp == 98


def test_blastn_from_outfmt_multiple_lines():
    data = "99.5,150,98\n88.0,120,90\n"
    results = Blastn.from_outfmt(data, BlastnOutFmt.FMT_1)
    assert len(results) == 2
    assert results[1].pident == 88.0
    assert results[1].length == 120


def test_blastn_from_outfmt_empty():
    assert Blastn.from_outfmt("", BlastnOutFmt.FMT_1) == []
