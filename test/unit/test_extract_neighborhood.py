"""Unit tests for the pure helpers in sarand.extract_neighborhood."""
import gfapy
import networkx as nx
import pytest

from sarand.extract_neighborhood import (
    calculate_coverage,
    create_directed_graph_nx,
    get_sequence_path,
    read_paths_file,
    reverse_complement,
    write_paths_file,
)


@pytest.mark.parametrize(
    "seq,expected",
    [
        ("AAAA", "TTTT"),
        ("ACGT", "ACGT"),
        ("ACCGGT", "ACCGGT"),
        ("GGGGCCCC", "GGGGCCCC"),
        ("ATGC", "GCAT"),
    ],
)
def test_reverse_complement(seq, expected):
    assert reverse_complement(seq) == expected


def test_reverse_complement_is_involution():
    seq = "ACGTTGCAACGT"
    assert reverse_complement(reverse_complement(seq)) == seq


class TestCalculateCoverage:
    def test_metaspades(self):
        node = {"KC": 80, "sequence": "AAAACCCC"}  # len 8, kmer 4 -> 80/4
        assert calculate_coverage(node, 4, "1+", "metaspades") == 20.0

    def test_bcalm_uses_same_formula(self):
        node = {"KC": 60, "sequence": "A" * 10}  # 60 / (10 - 4)
        assert calculate_coverage(node, 4, "1+", "bcalm") == 10.0

    def test_megahit_parses_from_node_name(self):
        assert calculate_coverage(None, 4, "k55_cov_12.5_x", "megahit") == 12.5

    def test_unknown_assembler_exits(self):
        with pytest.raises(SystemExit):
            calculate_coverage({"KC": 1, "sequence": "AA"}, 1, "n", "unknown")


class TestPathsFileRoundTrip:
    def test_write_then_read(self, tmp_path):
        paths = {
            ("1+", "2+"): "ACGTACGT",
            ("3-", "4+", "5-"): "TTTTGGGG",
        }
        f = tmp_path / "paths.txt"
        write_paths_file(f, paths)
        assert read_paths_file(f) == paths

    def test_append_mode(self, tmp_path):
        f = tmp_path / "paths.txt"
        write_paths_file(f, {("1+",): "AAAA"})
        write_paths_file(f, {("2+",): "CCCC"}, mode="a")
        assert read_paths_file(f) == {("1+",): "AAAA", ("2+",): "CCCC"}


class TestCreateDirectedGraph:
    def _graph(self):
        g = gfapy.Gfa()
        g.append("H\tVN:Z:1.0")
        g.append("S\t1\tAAAACCCC\tKC:i:80")
        g.append("S\t2\tCCCCGGGG\tKC:i:40")
        g.append("L\t1\t+\t2\t+\t4M")
        return create_directed_graph_nx(g)

    def test_creates_plus_and_minus_nodes(self):
        dg = self._graph()
        assert set(dg.nodes) == {"1+", "1-", "2+", "2-"}

    def test_minus_node_is_reverse_complement(self):
        dg = self._graph()
        assert dg.nodes["1+"]["sequence"] == "AAAACCCC"
        assert dg.nodes["1-"]["sequence"] == "GGGGTTTT"

    def test_forward_edge_overlap_and_weight(self):
        dg = self._graph()
        assert dg["1+"]["2+"]["overlap"] == 4
        # weight = len(dest seq) - overlap = 8 - 4
        assert dg["1+"]["2+"]["weight"] == 4

    def test_reverse_complement_edge_present(self):
        """A '1+ -> 2+' link implies a '2- -> 1-' link in the rc direction."""
        dg = self._graph()
        assert dg.has_edge("2-", "1-")


class TestGetSequencePath:
    def _graph(self):
        dg = nx.DiGraph()
        dg.add_node("a", sequence="AAAATTTT")
        dg.add_node("b", sequence="TTTTGGGG")
        dg.add_edge("a", "b", overlap=4)
        return dg

    def test_downstream_concatenation(self):
        dg = self._graph()
        # len_target keeps last 8 of node a (whole), then b minus the 4bp overlap
        seq = get_sequence_path(dg, ["a", "b"], threshold=100, len_target=8, up_down="down")
        assert seq == "AAAATTTT" + "GGGG"

    def test_downstream_threshold_trims(self):
        dg = self._graph()
        seq = get_sequence_path(dg, ["a", "b"], threshold=4, len_target=8, up_down="down")
        assert seq == "AAAA"

    def test_upstream_keeps_tail(self):
        dg = self._graph()
        seq = get_sequence_path(dg, ["a", "b"], threshold=100, len_target=8, up_down="up")
        # b (minus overlap on its tail) then a's leading len_target bases
        assert seq.endswith("AAAATTTT")
