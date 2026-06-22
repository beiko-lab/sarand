"""Unit tests for the FASTA I/O helpers in sarand.util.sequence.

The blastn-backed ``compare_two_sequences`` is exercised by the functional
test rather than here, since it shells out to the blastn binary.
"""
from sarand.util.sequence import create_fasta_file, retrieve_target


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
