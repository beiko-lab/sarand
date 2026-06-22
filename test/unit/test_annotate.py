"""Unit tests for sarand.util.annotate ORF calling and de-duplication.

``unnamed_genes_similar`` and the blastn-backed branches of
``annotations_identical`` shell out to blastn and are covered by the functional
test; here we exercise the named-gene and ORF-calling logic directly.
"""
from sarand.util.annotate import (
    annotations_identical,
    call_orfs,
    partition_genes_around_amr,
    similar_annotation_exists,
    write_orf_files,
)


def _gene(start, end, gene="", seq_value="", coverage=1.0):
    return {
        "gene": gene,
        "product": "",
        "length": abs(end - start) + 1,
        "start_pos": start,
        "end_pos": end,
        "coverage": coverage,
        "seq_value": seq_value,
        "seq_name": "s1",
        "target_amr": None,
    }


class TestCallOrfs:
    def test_finds_orf_and_populates_fields(self):
        # ATG + 40 codons + stop -> a single clean ORF
        seq = ("ATG" + "GCT" * 40 + "TAA").lower()
        result = call_orfs(seq)
        assert len(result) >= 1
        orf = result[0]
        for key in ("gene", "product", "length", "start_pos", "end_pos",
                    "strand", "nt_seq", "aa_seq", "seq_value"):
            assert key in orf
        # call_orfs uppercases internally; seq_value keeps the (rstripped) input
        assert orf["seq_value"] == seq

    def test_tolerates_trailing_newline(self):
        seq = "ATG" + "GCT" * 40 + "TAA" + "\n"
        result = call_orfs(seq)
        assert result[0]["seq_value"] == seq.rstrip("\n")

    def test_no_genes_returns_empty(self):
        assert call_orfs("ACAC") == []


class TestPartitionGenesAroundAmr:
    def test_partitions_upstream_target_downstream(self):
        # AMR is the lower-case middle region (idx 10..19)
        seq = "A" * 10 + "c" * 10 + "G" * 10
        seq_info = [
            _gene(1, 9, seq_value=seq),     # upstream
            _gene(11, 20, seq_value=seq),   # target (overlaps lower-case)
            _gene(21, 30, seq_value=seq),   # downstream
        ]
        found, amr, up, down, si = partition_genes_around_amr(seq, seq_info)
        assert found is True
        assert amr["start_pos"] == 11
        assert amr["target_amr"] == "yes"
        assert len(up) == 1 and up[0]["start_pos"] == 1
        assert len(down) == 1 and down[0]["start_pos"] == 21

    def test_amr_at_end_of_sequence(self):
        seq = "A" * 10 + "c" * 10  # no downstream upper-case region
        seq_info = [_gene(1, 9, seq_value=seq), _gene(11, 20, seq_value=seq)]
        found, amr, up, down, si = partition_genes_around_amr(seq, seq_info)
        assert found is True
        assert amr["start_pos"] == 11


class TestAnnotationDeduplication:
    def test_identical_named_annotations(self):
        a = [{"gene": "x", "start_pos": 1, "end_pos": 2},
             {"gene": "y", "start_pos": 3, "end_pos": 4}]
        b = [{"gene": "x", "start_pos": 1, "end_pos": 2},
             {"gene": "y", "start_pos": 3, "end_pos": 4}]
        assert annotations_identical(a, b, "/tmp") is True

    def test_different_gene_names_not_identical(self):
        a = [{"gene": "x", "start_pos": 1, "end_pos": 2}]
        b = [{"gene": "z", "start_pos": 1, "end_pos": 2}]
        assert annotations_identical(a, b, "/tmp") is False

    def test_different_lengths_not_identical(self):
        a = [{"gene": "x", "start_pos": 1, "end_pos": 2}]
        b = [{"gene": "x", "start_pos": 1, "end_pos": 2},
             {"gene": "y", "start_pos": 3, "end_pos": 4}]
        assert annotations_identical(a, b, "/tmp") is False

    def test_similar_annotation_exists(self):
        new = [{"gene": "x", "start_pos": 1, "end_pos": 2}]
        existing = [
            [{"gene": "q", "start_pos": 1, "end_pos": 2}],
            [{"gene": "x", "start_pos": 1, "end_pos": 2}],
        ]
        assert similar_annotation_exists(new, existing, "/tmp") is True

    def test_similar_annotation_absent(self):
        new = [{"gene": "x", "start_pos": 1, "end_pos": 2}]
        existing = [[{"gene": "q", "start_pos": 1, "end_pos": 2}]]
        assert similar_annotation_exists(new, existing, "/tmp") is False


def test_write_orf_files(tmp_path):
    seq = "ATG" + "GCT" * 40 + "TAA"
    orfs = call_orfs(seq)
    for orf in orfs:
        orf["seq_name"] = "neighborhood_1"
    ffn, faa, gff = write_orf_files([orfs], tmp_path / "orfs")
    assert ffn.suffix == ".ffn" and faa.suffix == ".faa" and gff.suffix == ".gff"
    assert gff.read_text().startswith("##gff-version 3\n")
    # one CDS line per ORF plus the header line
    cds_lines = [l for l in gff.read_text().splitlines() if "\tCDS\t" in l]
    assert len(cds_lines) == len(orfs)
    assert faa.read_text().startswith(">neighborhood_1_orf1")
