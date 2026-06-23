"""Helpers for deriving and normalising target sequence and file names."""
from __future__ import annotations

from pathlib import Path


def extract_name_from_file_name(file_name: str | Path) -> str:
    """Return the base name of a file without its directory or extension."""
    return Path(file_name).stem


def target_name_from_comment(target_comment: str) -> str:
    """Derive a normalised sequence name from a FASTA comment/header line.

    Only the portion of the header up to the first whitespace is considered,
    mirroring the way BLAST (and therefore Bandage) truncates query names at the
    first whitespace. This keeps the name derived here consistent with the target
    name parsed back out of the Bandage alignment output. Without it, databases
    whose headers carry a trailing token after the gene name -- e.g. the NCBI
    ``AMR_CDS.fa`` headers ``...|geneName accession:coords`` -- would have the
    trailing ``accession:coords`` folded into the name here but stripped by
    Bandage, so the two names never match and every hit is silently discarded.
    """
    header = target_comment.split()
    prefix = header[0] if header else ""
    return (
        prefix.split("|")[-1]
        .strip()
        .replace(" ", "_")
        .replace("'", ";")
        .replace("/", "]")
    )


def restricted_target_name_from_modified_name(target_name: str) -> str:
    """Convert a target name into a filesystem-safe form."""
    restricted = target_name.replace(";", "SS")
    return "".join(
        e for e in restricted if e.isalpha() or e.isnumeric() or e == "_" or e == "-"
    )
