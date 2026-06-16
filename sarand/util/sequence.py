"""FASTA I/O helpers and blastn-based sequence comparison."""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict

from Bio import SeqIO

from sarand.external.blastn import Blastn
from sarand.model.fasta_seq import FastaSeq
from sarand.util.naming import target_name_from_comment


def create_fasta_file(seq: str, output_dir: str | Path, comment: str = "> sequence:\n",
                      file_name: str = "temp") -> str:
    """
    To create a fasta file for a sequence
    Parameters:
        seq: the sequence to be written into the file
        output_dir: the output directory address
        comment: the comment to be written into fasta file
        file_name: the name of the fasta file
    Return:
        the address of the fasta file
    """
    myfile_name = Path(output_dir) / (file_name + ".fasta")
    if myfile_name.is_file():
        myfile_name.unlink()
    with open(myfile_name, 'w') as myfile:
        myfile.write(comment)
        if not comment.endswith("\n"):
            myfile.write("\n")
        myfile.write(seq)
        if not seq.endswith("\n"):
            myfile.write("\n")
    # Returned as a plain string: callers concatenate this path into messages
    # and pass it to subprocesses.
    return str(myfile_name)


def retrieve_target(file_path: str | Path) -> tuple[str, str]:
    """
    To read the target gene from the text file.
    Parameters:
        file_path:	the address of the file containing the AMR gene
    Return:
        the sequence of the AMR gene in lower case
    """
    target_name = ""
    with open(file_path) as fp:
        for line in fp:
            # skip comment line
            if line.startswith(">"):
                target_name = target_name_from_comment(line[:-1])
                continue
            return line, target_name


def compare_two_sequences(
        subject: str,
        query: str,
        output_dir: str | Path,
        threshold: int = 90,
        switch_allowed: bool = True,
        return_file: bool = False,
        subject_coverage: bool = True,
        blast_ext: str = "",
) -> bool:
    """
    To compare one sequence (shorter sequence) against the other one (longer
    sequence) using blastn.
    """
    # make sure subject is the longer sequence
    if switch_allowed and len(subject) < len(query):
        subject, query = query, subject

    # Write to a temporary directory to prevent any race condition when this is
    # called via the recursion of the neighbourhood extraction.
    with tempfile.TemporaryDirectory() as tmp_dir:
        query_file_name = Path(tmp_dir) / "query.fasta"
        with open(query_file_name, "w") as query_file:
            query_file.write("> query \n")
            query_file.write(query)
        subject_file_name = Path(tmp_dir) / "subject.fasta"
        with open(subject_file_name, "w") as subject_file:
            subject_file.write("> subject \n")
            subject_file.write(subject)

        blastn = Blastn.run_for_sarand_compare_two_sequences(
            query=query_file_name,
            subject=subject_file_name
        )

        if return_file:
            raise NotImplementedError('Unable to return file.')

        for row in blastn.results:
            identity = int(row.pident)
            coverage = int(row.length / len(subject) * 100)
            q_coverage = row.qcovhsp

            if subject_coverage and identity >= threshold and coverage >= threshold:
                return True
            if not subject_coverage and identity >= threshold and q_coverage >= threshold:
                return True
    return False


def extract_amr_sequences(path: Path) -> Dict[str, FastaSeq]:
    """Extract the AMR sequences from a FASTA file."""
    out = dict()
    with path.open() as f:
        for record in SeqIO.parse(f, "fasta"):
            target_name = target_name_from_comment(record.description)
            if target_name in out:
                raise ValueError(f"Duplicate AMR name {target_name} in {path}")
            out[target_name] = FastaSeq(
                seq=str(record.seq),
                fasta_id=str(record.description)
            )
    return out
