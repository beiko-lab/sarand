"""Unit tests for sarand.util.graph_path grouping/filtering of Bandage results."""
from sarand.external.bandage import BandageResult
from sarand.util.graph_path import (
    read_path_info_from_align_file_with_multiple_targets,
)


def make_result(query, path="(10) 1+, 2+ (20)", covered="99%", identity="99%"):
    return BandageResult(
        query=query,
        path_with_start_end=path,
        length=100,
        query_covered_by_path=covered,
        query_covered_by_hits=covered,
        mean_hit_identity=identity,
        total_hit_mismatches=0,
        total_hit_gap_opens=0,
        relative_length="100%",
        length_discrepancy=0,
        e_value_product=1e-50,
        sequence="ACGT",
    )


def test_groups_passing_results_by_target(tmp_path):
    results = [
        make_result("gb|X|ARO|ErmA"),
        make_result("gb|X|ARO|ErmA", path="(5) 3+ (15)"),
        make_result("gb|X|ARO|ANT9-Ia"),
    ]
    grouped = read_path_info_from_align_file_with_multiple_targets(
        tmp_path, results, min_target_identity=95, min_target_coverage=95
    )
    assert set(grouped) == {"ErmA", "ANT9-Ia"}
    assert len(grouped["ErmA"]) == 2
    first = grouped["ErmA"][0]
    assert first["nodes"] == ["1", "2"]
    assert first["orientations"] == ["+", "+"]
    assert first["start_pos"] == 10
    assert first["end_pos"] == 20


def test_filters_out_below_threshold(tmp_path):
    results = [
        make_result("gb|X|ARO|ErmA", covered="50%", identity="99%"),  # low coverage
        make_result("gb|X|ARO|ErmA", covered="99%", identity="50%"),  # low identity
    ]
    grouped = read_path_info_from_align_file_with_multiple_targets(
        tmp_path, results, min_target_identity=95, min_target_coverage=95
    )
    assert grouped == {}


def test_debug_writes_json(tmp_path):
    results = [make_result("gb|X|ARO|ErmA")]
    read_path_info_from_align_file_with_multiple_targets(
        tmp_path, results, debug=True
    )
    assert (tmp_path / "debug_alignment_coverage_identity.json").is_file()
