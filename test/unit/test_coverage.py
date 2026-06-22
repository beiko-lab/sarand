"""Unit tests for sarand.coverage pure calculation/formatting helpers."""
import csv

from sarand.coverage import (
    _format_coverage,
    _gene_path,
    _target_coverage,
    _target_gene_coords,
    find_amr_coverage_in_seq,
    find_gene_coverage,
    read_path_coverage_info,
)


class TestReadPathCoverageInfo:
    def _write(self, tmp_path, rows):
        f = tmp_path / "paths_info.csv"
        with open(f, "w", newline="") as fh:
            w = csv.DictWriter(
                fh, fieldnames=["sequence", "node", "coverage", "start", "end"]
            )
            w.writeheader()
            for r in rows:
                w.writerow(r)
        return f

    def test_groups_nodes_by_sequence(self, tmp_path):
        f = self._write(
            tmp_path,
            [
                {"sequence": "s1", "node": "1", "coverage": "10", "start": "0", "end": "5"},
                {"sequence": "s1", "node": "2", "coverage": "12", "start": "6", "end": "9"},
                {"sequence": "s2", "node": "3", "coverage": "20", "start": "0", "end": "4"},
            ],
        )
        result = read_path_coverage_info(f)
        assert len(result) == 2
        assert len(result[0]) == 2
        assert result[0][0] == {"node": "1", "coverage": 10.0, "start": 0, "end": 5}
        assert result[1][0]["node"] == "3"

    def test_skips_duplicate_header_rows(self, tmp_path):
        f = self._write(
            tmp_path,
            [
                {"sequence": "s1", "node": "1", "coverage": "10", "start": "0", "end": "5"},
                {"sequence": "coverage", "node": "node", "coverage": "coverage",
                 "start": "start", "end": "end"},
            ],
        )
        result = read_path_coverage_info(f)
        assert len(result) == 1
        assert len(result[0]) == 1


def test_find_gene_coverage_single_node():
    # gene spans 1-based positions 2..5 -> 0-based 1..4 (4 bases)
    seq_info_list = [{"start_pos": 2, "end_pos": 5}]
    path_info = [{"start": 0, "end": 10, "coverage": 10.0}]
    assert find_gene_coverage(seq_info_list, path_info) == [10.0]


def test_find_gene_coverage_spans_two_nodes():
    # gene 0-based 1..6 spanning node A (0..3 cov 10) and node B (4..10 cov 20)
    seq_info_list = [{"start_pos": 2, "end_pos": 7}]
    path_info = [
        {"start": 0, "end": 3, "coverage": 10.0},
        {"start": 4, "end": 10, "coverage": 20.0},
    ]
    # bases 1..3 (3 bases) at cov 10 + bases 4..6 (3 bases) at cov 20
    # = (3*10 + 3*20) / 6 = 90/6 = 15
    assert find_gene_coverage(seq_info_list, path_info) == [15.0]


def test_find_amr_coverage_in_seq():
    seq = "A" * 10 + "c" * 10 + "G" * 10  # AMR is lower-case region (idx 10..19)
    seq_info = [
        {"start_pos": 1, "end_pos": 9, "seq_value": seq, "coverage": 5.0},
        {"start_pos": 11, "end_pos": 20, "seq_value": seq, "coverage": 42.0},
        {"start_pos": 21, "end_pos": 30, "seq_value": seq, "coverage": 7.0},
    ]
    coverage, index, error = find_amr_coverage_in_seq(seq_info)
    assert error is False
    assert index == 1
    assert coverage == 42.0


def test_find_amr_coverage_no_overlap_returns_error():
    seq = "A" * 10 + "c" * 10 + "G" * 10
    seq_info = [
        {"start_pos": 1, "end_pos": 5, "seq_value": seq, "coverage": 5.0},
    ]
    coverage, index, error = find_amr_coverage_in_seq(seq_info)
    assert error is True
    assert index == -1


def test_gene_path_orders_lo_hi():
    seq_info = [
        {"start_pos": 5, "end_pos": 1},  # reversed coords normalised to lo-hi
        {"start_pos": 10, "end_pos": 20},
    ]
    assert _gene_path(seq_info) == "1-5;10-20"


def test_format_coverage():
    assert _format_coverage(None) == ""
    assert _format_coverage("") == ""
    assert _format_coverage(12.345) == "12.35"
    assert _format_coverage("9.0") == "9.0"


def test_target_coverage_and_coords():
    seq_info = [
        {"start_pos": 1, "end_pos": 9, "coverage": 5.0, "target_amr": None},
        {"start_pos": 20, "end_pos": 11, "coverage": 42.0, "target_amr": "yes"},
    ]
    assert _target_coverage(seq_info) == 42.0
    assert _target_gene_coords(seq_info) == "11-20"


def test_target_coverage_missing_returns_blank():
    seq_info = [{"start_pos": 1, "end_pos": 9, "coverage": 5.0, "target_amr": None}]
    assert _target_coverage(seq_info) == ""
    assert _target_gene_coords(seq_info) == ""
