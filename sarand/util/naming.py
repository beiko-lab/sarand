"""Helpers for deriving and normalising AMR and file names."""
import os


def extract_name_from_file_name(file_name):
    """Return the base name of a file without its directory or extension."""
    return os.path.splitext(os.path.basename(file_name))[0]


def amr_name_from_comment(amr_comment):
    """Derive a normalised AMR name from a FASTA comment/header line."""
    return (
        amr_comment.split("[")[0]
        .split("|")[-1]
        .strip()
        .replace(" ", "_")
        .replace("'", ";")
        .replace("/", "]")
    )


def restricted_amr_name_from_modified_name(amr_name):
    """Convert an AMR name into a filesystem-safe form."""
    restricted = amr_name.replace(";", "SS")
    return "".join(
        e for e in restricted if e.isalpha() or e.isnumeric() or e == "_" or e == "-"
    )
