"""Unit tests for sarand.external.bandage parsing and command building."""
from pathlib import Path

import pytest

from sarand.external.bandage import Bandage, BandageParams, BandageResult


def make_result(query="q", path="(1363) 69-, 2193+ (1786)",
                covered_by_path="90%", mean_identity="99.5%"):
    return BandageResult(
        query=query,
        path_with_start_end=path,
        length=100,
        query_covered_by_path=covered_by_path,
        query_covered_by_hits="90%",
        mean_hit_identity=mean_identity,
        total_hit_mismatches=0,
        total_hit_gap_opens=0,
        relative_length="100%",
        length_discrepancy=0,
        e_value_product=1e-50,
        sequence="ACGT",
    )


class TestBandageResult:
    def test_path_strips_position_annotations(self):
        assert make_result().path == "69-, 2193+"

    def test_path_start_and_end(self):
        r = make_result()
        assert r.path_start == 1363
        assert r.path_end == 1786

    def test_path_start_end_default_zero_when_unannotated(self):
        r = make_result(path="69-, 2193+")
        assert r.path_start == 0
        assert r.path_end == 0

    def test_path_to_sarand_nodes_and_orientations(self):
        nodes, orientations = make_result().path_to_sarand
        assert nodes == ["69", "2193"]
        assert orientations == ["-", "+"]

    def test_path_to_sarand_empty_path_raises(self):
        with pytest.raises(Exception):
            make_result(path="()").path_to_sarand

    def test_coverage_and_identity_pct(self):
        r = make_result(covered_by_path="88.5%", mean_identity="97.2%")
        assert r.coverage_pct == 88.5
        assert r.identity_pct == 97.2

    def test_target_name_normalises_query(self):
        r = make_result(query="gb|X|ARO|aac(6')-Ie")
        assert r.target_name == "aac6SS-Ie"

    def test_repr_is_query(self):
        assert repr(make_result(query="myquery")) == "myquery"


class TestBandageParams:
    def test_as_cmd_minimal(self):
        params = BandageParams(
            graph=Path("g.gfa"),
            reads=Path("r.fasta"),
            outputfile=Path("out.tsv"),
        )
        cmd = [str(x) for x in params.as_cmd()]
        assert cmd[:2] == ["Bandage", "querypaths"]
        assert str(Path("g.gfa").absolute()) in cmd
        assert str(Path("r.fasta").absolute()) in cmd

    def test_as_cmd_includes_optional_params(self):
        params = BandageParams(
            graph=Path("g.gfa"),
            reads=Path("r.fasta"),
            outputfile=Path("out.tsv"),
            pathnodes=6,
            minpatcov=0.9,
            minmeanid=0.5,
            minhitcov=0.9,
            verbose=True,
        )
        cmd = [str(x) for x in params.as_cmd()]
        assert "--pathnodes" in cmd and "6" in cmd
        assert "--minpatcov" in cmd and "0.9" in cmd
        assert "--minmeanid" in cmd and "0.5" in cmd
        assert "--verbose" in cmd

    def test_from_cli_args_pairs_and_flags(self):
        params = BandageParams.from_cli_args(
            [["pathnodes", "10"], ["minmeanid", "0.7"], ["verbose"]]
        )
        assert params.pathnodes == 10
        assert params.minmeanid == 0.7
        assert params.verbose is True

    def test_from_cli_args_none(self):
        params = BandageParams.from_cli_args(None)
        assert params.pathnodes is None

    def test_update_from_dictionary_casts_types(self):
        params = BandageParams()
        params.update_from_dictionary(
            {"pathnodes": "8", "minpatcov": "0.85", "minlendis": "-5"}
        )
        assert params.pathnodes == 8 and isinstance(params.pathnodes, int)
        assert params.minpatcov == 0.85 and isinstance(params.minpatcov, float)
        assert params.minlendis == -5

    def test_update_from_object_takes_non_none(self):
        base = BandageParams(pathnodes=6, minmeanid=0.5)
        other = BandageParams(pathnodes=10)  # minmeanid stays None
        base.update_from_object(other)
        assert base.pathnodes == 10
        assert base.minmeanid == 0.5  # unchanged because other's was None


class TestBandageReadFile:
    def test_read_file_parses_rows(self, tmp_path):
        tsv = tmp_path / "bandage.tsv"
        header = "\t".join(
            [
                "Query", "Path", "Length", "QueryCovByPath", "QueryCovByHits",
                "MeanHitId", "Mismatches", "GapOpens", "RelLength",
                "LengthDisc", "EValue", "Sequence",
            ]
        )
        row = "\t".join(
            [
                "aac(6')-Ie", "(1363) 69-, 2193+ (1786)", "100", "90%", "90%",
                "99.5%", "0", "0", "100%", "0", "1e-50", "ACGT",
            ]
        )
        tsv.write_text(header + "\n" + row + "\n")
        results = Bandage.read_file(tsv)
        assert len(results) == 1
        assert results[0].length == 100
        assert results[0].coverage_pct == 90.0
        assert results[0].path_to_sarand == (["69", "2193"], ["-", "+"])

    def test_read_file_header_only_returns_empty(self, tmp_path):
        tsv = tmp_path / "bandage.tsv"
        tsv.write_text("Query\tPath\n")
        assert Bandage.read_file(tsv) == []
