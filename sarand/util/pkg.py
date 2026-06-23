import sys
from pathlib import Path

# Filenames of the reference FASTAs bundled under ``sarand/data``.
_BUNDLED_FASTA = {
    "card": "CARD_AMR_seq.fasta",
    "ncbi": "NCBI_AMR_seq.fasta",
}


def get_pkg_data_dir() -> Path:
    """Return the package ``data`` directory."""
    return Path(sys.modules['sarand'].__file__).parent / 'data'


def get_pkg_fasta_path(database: str = "card") -> Path:
    """Return the path to the bundled reference FASTA for ``database``.

    https://docs.python.org/3/library/pkgutil.html#pkgutil.get_data
    """
    try:
        name = _BUNDLED_FASTA[database]
    except KeyError:
        raise ValueError(f"No bundled FASTA for database {database!r}")
    return get_pkg_data_dir() / name


def get_pkg_card_fasta_path() -> Path:
    """Return the path to the bundled CARD FASTA (backwards-compatible alias)."""
    return get_pkg_fasta_path("card")
