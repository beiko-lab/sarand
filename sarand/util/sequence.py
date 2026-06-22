"""FASTA I/O helpers and minimap2-based sequence comparison."""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict, Set, Tuple

from sarand.external.minimap2 import Minimap2
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
    Read the target gene sequence and name from a FASTA file.
    Parameters:
        file_path:	the address of the file containing the target gene
    Return:
        (sequence, target_name): the first sequence line and the name parsed from
        the FASTA header
    """
    target_name = ""
    with open(file_path) as fp:
        for line in fp:
            # skip comment line
            if line.startswith(">"):
                target_name = target_name_from_comment(line[:-1])
                continue
            return line, target_name


def similar_sequence_pairs(
        queries: Dict[str, str],
        subjects: Dict[str, str],
        threshold: int = 90,
) -> Set[Tuple[str, str]]:
    """
    Find which query sequences are near-identical to which subject sequences,
    using a single batched minimap2 alignment of all queries against all
    subjects (rather than one alignment per pair).

    A (query, subject) pair is considered similar when an alignment between them
    covers at least ``threshold`` percent of the longer sequence at at least
    ``threshold`` percent identity.
    Parameters:
        queries: mapping of query id -> nucleotide sequence
        subjects: mapping of subject id -> nucleotide sequence
        threshold: minimum percent identity and coverage to call a pair similar
    Return:
        the set of (query_id, subject_id) pairs that are similar
    """
    # drop empty sequences (nothing to align) and short-circuit empty inputs
    queries = {k: v for k, v in queries.items() if v}
    subjects = {k: v for k, v in subjects.items() if v}
    if not queries or not subjects:
        return set()

    # Write to a temporary directory to avoid any race condition when this is
    # called from the parallel per-target neighborhood processing.
    with tempfile.TemporaryDirectory() as tmp_dir:
        query_file = Path(tmp_dir) / "queries.fasta"
        with open(query_file, "w") as fh:
            for qid, seq in queries.items():
                fh.write(f">{qid}\n{seq}\n")
        subject_file = Path(tmp_dir) / "subjects.fasta"
        with open(subject_file, "w") as fh:
            for sid, seq in subjects.items():
                fh.write(f">{sid}\n{seq}\n")

        results = Minimap2.run(query=query_file, target=subject_file)

    matched: Set[Tuple[str, str]] = set()
    for row in results:
        if row.identity_pct >= threshold and row.coverage_pct >= threshold:
            matched.add((row.query, row.target))
    return matched
