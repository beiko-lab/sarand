class FastaSeq:
    """Wrapper for a sequence in a fasta file."""

    __slots__ = ('seq', 'fasta_id')

    def __init__(self, seq: str, fasta_id: str):
        """
        Parameters:
            seq: The sequence.
            fasta_id: The id in the fasta file.
        """
        self.seq: str = seq
        self.fasta_id: str = fasta_id
