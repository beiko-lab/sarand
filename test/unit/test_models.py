"""Unit tests for the FastaSeq model."""
import pytest

from sarand.model.fasta_seq import FastaSeq


def test_fasta_seq_stores_fields():
    fs = FastaSeq(seq="ACGT", fasta_id="seq1")
    assert fs.seq == "ACGT"
    assert fs.fasta_id == "seq1"


def test_fasta_seq_is_slotted():
    """__slots__ should prevent arbitrary attribute assignment."""
    fs = FastaSeq("ACGT", "seq1")
    assert not hasattr(fs, "__dict__")
    with pytest.raises(AttributeError):
        fs.unexpected = "value"
