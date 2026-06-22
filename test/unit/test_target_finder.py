"""Unit tests for the path-overlap logic in sarand.target_finder."""
from sarand.target_finder import target_path_overlap


def _path(nodes, orientations, start_pos, end_pos):
    return {
        "nodes": nodes,
        "orientations": orientations,
        "start_pos": start_pos,
        "end_pos": end_pos,
    }


def test_all_new_paths_overlap_existing():
    found = [[_path(["1"], ["+"], 0, 100)]]
    new = [_path(["1"], ["+"], 0, 100)]
    all_overlap, id_list = target_path_overlap(found, new, new_target_len=100)
    assert all_overlap is True
    assert id_list == [0]


def test_unique_new_path_when_nodes_differ():
    found = [[_path(["1"], ["+"], 0, 100)]]
    new = [_path(["2"], ["+"], 0, 100)]  # different node, no overlap
    all_overlap, id_list = target_path_overlap(found, new, new_target_len=100)
    assert all_overlap is False
    assert id_list is None


def test_partial_overlap_below_threshold_is_unique():
    found = [[_path(["1"], ["+"], 50, 60)]]
    # large coordinate discrepancy -> percent below the 95 default
    new = [_path(["1"], ["+"], 0, 100)]
    all_overlap, id_list = target_path_overlap(found, new, new_target_len=100)
    assert all_overlap is False


def test_orientation_must_match():
    found = [[_path(["1"], ["+"], 0, 100)]]
    new = [_path(["1"], ["-"], 0, 100)]
    all_overlap, id_list = target_path_overlap(found, new, new_target_len=100)
    assert all_overlap is False


def test_multiple_new_paths_all_must_overlap():
    found = [[_path(["1"], ["+"], 0, 100)], [_path(["2"], ["+"], 0, 100)]]
    new = [_path(["1"], ["+"], 0, 100), _path(["3"], ["+"], 0, 100)]
    # second new path (node 3) has no match -> not all overlap
    all_overlap, id_list = target_path_overlap(found, new, new_target_len=100)
    assert all_overlap is False
