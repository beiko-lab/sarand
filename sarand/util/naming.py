"""Helpers for deriving and normalising AMR and file names."""
from __future__ import annotations

from pathlib import Path


def extract_name_from_file_name(file_name: str | Path) -> str:
    """Return the base name of a file without its directory or extension."""
    return Path(file_name).stem


def amr_name_from_comment(amr_comment: str) -> str:
    """Derive a normalised AMR name from a FASTA comment/header line."""
    return (
        amr_comment.split("[")[0]
        .split("|")[-1]
        .strip()
        .replace(" ", "_")
        .replace("'", ";")
        .replace("/", "]")
    )


def restricted_amr_name_from_modified_name(amr_name: str) -> str:
    """Convert an AMR name into a filesystem-safe form."""
    restricted = amr_name.replace(";", "SS")
    return "".join(
        e for e in restricted if e.isalpha() or e.isnumeric() or e == "_" or e == "-"
    )
