"""Unit tests for the FASTA I/O helpers in sarand.util.sequence.

``similar_sequence_pairs`` shells out to minimap2; the tests that exercise it are
skipped when the binary is unavailable, and the path is also covered end-to-end
by the functional test.
"""
import random
import shutil

import pytest

from sarand.util.sequence import (
    create_fasta_file,
    retrieve_target,
    similar_sequence_pairs,
)

requires_minimap2 = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 not installed"
)


def _random_seq(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _with_substitutions(seq, positions):
    flip = {"A": "C", "C": "A", "G": "T", "T": "G"}
    chars = list(seq)
    for i in positions:
        chars[i] = flip[chars[i]]
    return "".join(chars)


def test_create_fasta_file_writes_comment_and_seq(tmp_path):
    path = create_fasta_file("ACGT", tmp_path, comment="> mygene\n", file_name="g")
    assert path == str(tmp_path / "g.fasta")
    assert open(path).read() == "> mygene\nACGT\n"


def test_create_fasta_file_adds_missing_newlines(tmp_path):
    """A comment/sequence without trailing newlines is normalised."""
    path = create_fasta_file("ACGT", tmp_path, comment="> mygene", file_name="g")
    assert open(path).read() == "> mygene\nACGT\n"


def test_create_fasta_file_overwrites_existing(tmp_path):
    create_fasta_file("AAAA", tmp_path, comment="> a\n", file_name="g")
    path = create_fasta_file("CCCC", tmp_path, comment="> c\n", file_name="g")
    assert open(path).read() == "> c\nCCCC\n"


def test_retrieve_target_parses_name_and_sequence(tmp_path):
    f = tmp_path / "target.fasta"
    f.write_text(">gb|X|ARO|aac(6')-Ie [Staph]\nACGTACGT\n")
    seq, name = retrieve_target(f)
    assert seq == "ACGTACGT\n"
    assert name == "aac(6;)-Ie"


class TestSimilarSequencePairs:
    def test_empty_inputs_return_empty(self):
        # no aligner is invoked when either side is empty
        assert similar_sequence_pairs({"q": "ACGT"}, {}) == set()
        assert similar_sequence_pairs({}, {"s": "ACGT"}) == set()

    @requires_minimap2
    def test_batched_match_and_mismatch(self):
        base = _random_seq(400, seed=7)
        near = _with_substitutions(base, [50, 200, 350])  # ~99% identity
        unrelated = _random_seq(400, seed=99)
        queries = {"q_match": near, "q_diff": unrelated}
        subjects = {"s_base": base}
        matched = similar_sequence_pairs(queries, subjects, threshold=90)
        assert ("q_match", "s_base") in matched
        assert ("q_diff", "s_base") not in matched

    @requires_minimap2
    def test_one_call_resolves_many_pairs(self):
        a = _random_seq(300, seed=1)
        b = _random_seq(300, seed=2)
        queries = {"qa": _with_substitutions(a, [10]), "qb": _with_substitutions(b, [20])}
        subjects = {"sa": a, "sb": b}
        matched = similar_sequence_pairs(queries, subjects, threshold=90)
        assert ("qa", "sa") in matched
        assert ("qb", "sb") in matched
        assert ("qa", "sb") not in matched
        assert ("qb", "sa") not in matched
