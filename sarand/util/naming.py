"""Helpers for deriving and normalising target sequence and file names."""
from __future__ import annotations

from pathlib import Path


def extract_name_from_file_name(file_name: str | Path) -> str:
    """Return the base name of a file without its directory or extension."""
    return Path(file_name).stem


def target_name_from_comment(target_comment: str) -> str:
    """Derive a normalised sequence name from a FASTA comment/header line."""
    return (
        target_comment.split("[")[0]
        .split("|")[-1]
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
