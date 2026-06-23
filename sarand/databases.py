"""Download, update and resolve the target-gene reference databases.

Two reference databases of nucleotide target genes are supported:

- ``card``  : the CARD "protein homolog model" nucleotide sequences
              (https://card.mcmaster.ca). This is the default.
- ``ncbi``  : the NCBI AMRFinderPlus reference gene CDS (``AMR_CDS.fa``)
              (https://ftp.ncbi.nlm.nih.gov/pathogen/...).

A copy of both databases is bundled with the package, so neither requires a
download before first use. ``sarand --update`` refreshes the cached copies to
the latest upstream releases.

Downloaded databases are cached in a user-writable directory (see
``get_db_dir``) and an updated copy is preferred over the bundled one. The
graph search aligns nucleotide target genes with Bandage+BLAST, so the
nucleotide files (CARD ``nucleotide_fasta_protein_homolog_model.fasta`` /
NCBI ``AMR_CDS.fa``) are used. They are normalised to one sequence per line on
install because the downstream parser reads the target FASTA line-by-line.
"""
from __future__ import annotations

import io
import json
import os
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Optional, TextIO

from Bio import SeqIO

from sarand.util.logger import LOG
from sarand.util.pkg import get_pkg_fasta_path, get_pkg_version

DATABASES = ("card", "ncbi")

CARD_DATA_URL = "https://card.mcmaster.ca/latest/data"
CARD_FASTA_MEMBER = "nucleotide_fasta_protein_homolog_model.fasta"

NCBI_BASE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/"
    "AMRFinderPlus/database/latest/"
)
NCBI_CDS_FILE = "AMR_CDS.fa"

_FASTA_NAME = {"card": "card_protein_homolog.fasta", "ncbi": "AMR_CDS.fasta"}


def get_db_dir() -> Path:
    """Return the user-writable directory used to cache downloaded databases.

    Overridable with ``$SARAND_DB_DIR``; otherwise ``$XDG_DATA_HOME/sarand`` or
    ``~/.local/share/sarand``.
    """
    override = os.environ.get("SARAND_DB_DIR")
    if override:
        return Path(override)
    xdg = os.environ.get("XDG_DATA_HOME")
    root = Path(xdg) if xdg else Path.home() / ".local" / "share"
    return root / "sarand"


def _db_subdir(database: str) -> Path:
    """Return the cache subdirectory for ``database``."""
    return get_db_dir() / database


def _fasta_path(database: str) -> Path:
    """Return the path of the installed FASTA for ``database``."""
    return _db_subdir(database) / _FASTA_NAME[database]


def _meta_path(database: str) -> Path:
    """Return the path of the version/metadata file for ``database``."""
    return _db_subdir(database) / "metadata.json"


def _load_meta(database: str) -> dict:
    """Read the cached metadata dict for ``database`` (empty dict if absent)."""
    path = _meta_path(database)
    if path.is_file():
        try:
            return json.loads(path.read_text())
        except (ValueError, OSError):
            return {}
    return {}


def _save_meta(database: str, meta: dict) -> None:
    """Write the metadata dict for ``database`` to its cache directory."""
    _meta_path(database).write_text(json.dumps(meta, indent=2) + "\n")


def get_target_fasta(database: str = "card") -> Path:
    """Resolve the target-gene FASTA for ``database``.

    A downloaded/updated copy is preferred; otherwise we fall back to the copy
    bundled with the package (both ``card`` and ``ncbi`` are bundled). The
    resolved database, version and source are logged.
    """
    if database not in DATABASES:
        raise ValueError(f"Unknown database {database!r}; choose from {DATABASES}")
    fasta = _fasta_path(database)
    if fasta.is_file():
        version = _load_meta(database).get("version", "unknown")
        LOG.info(f"Using {database.upper()} database version {version} (updated copy: {fasta})")
        return fasta
    bundled = get_pkg_fasta_path(database)
    if bundled.is_file():
        version = get_pkg_version(database)
        LOG.info(f"Using {database.upper()} database version {version} (bundled)")
        return bundled
    raise FileNotFoundError(
        f"No local '{database}' database found. "
        f"Download it first with: sarand --update --database {database}"
    )


def _normalise_fasta(src_handle: TextIO, dest: Path) -> int:
    """Write records from ``src_handle`` to ``dest`` as one sequence per line.

    The downstream target-gene parser reads the FASTA line-by-line (one
    sequence per record line), so any wrapped input must be unwrapped here.
    Written to a temporary file first and moved into place atomically.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    count = 0
    with tmp.open("w") as out:
        for record in SeqIO.parse(src_handle, "fasta"):
            out.write(f">{record.description}\n{record.seq}\n")
            count += 1
    tmp.replace(dest)
    return count


def update_card(force: bool = False) -> Path:
    """Download/refresh the CARD protein-homolog nucleotide database if newer."""
    fasta = _fasta_path("card")
    meta = _load_meta("card")

    LOG.info("Checking CARD for updates...")
    head = urllib.request.Request(CARD_DATA_URL, method="HEAD")
    with urllib.request.urlopen(head, timeout=60) as resp:
        last_modified = resp.headers.get("Last-Modified")

    if (
        not force
        and fasta.is_file()
        and last_modified is not None
        and meta.get("last_modified") == last_modified
    ):
        LOG.info(f"CARD is already up to date (v{meta.get('version', 'unknown')}).")
        return fasta

    LOG.info("Downloading CARD data...")
    with urllib.request.urlopen(CARD_DATA_URL, timeout=300) as resp:
        archive = resp.read()
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:bz2") as tf:
        version = "unknown"
        card_json = _find_member(tf, "card.json")
        if card_json is not None:
            try:
                version = str(json.load(tf.extractfile(card_json)).get("_version", "unknown"))
            except (ValueError, OSError):
                pass
        member = _find_member(tf, CARD_FASTA_MEMBER)
        if member is None:
            raise RuntimeError(f"{CARD_FASTA_MEMBER} not found in CARD archive")
        LOG.info(f"Installing CARD v{version} (previously: {meta.get('version', 'none')})")
        n = _normalise_fasta(io.TextIOWrapper(tf.extractfile(member)), fasta)

    _save_meta("card", {"version": version, "last_modified": last_modified})
    LOG.info(f"Installed {n} CARD target genes to {fasta}")
    return fasta


def update_ncbi(force: bool = False) -> Path:
    """Download/refresh the NCBI AMRFinderPlus CDS database if newer."""
    fasta = _fasta_path("ncbi")
    meta = _load_meta("ncbi")

    LOG.info("Checking NCBI AMRFinderPlus for updates...")
    with urllib.request.urlopen(NCBI_BASE_URL + "version.txt", timeout=60) as resp:
        version = resp.read().decode().strip()

    if not force and fasta.is_file() and meta.get("version") == version:
        LOG.info(f"NCBI database is already up to date ({version}).")
        return fasta

    LOG.info(f"Installing NCBI AMRFinderPlus {version} (previously: {meta.get('version', 'none')})")
    # Download to a temp file first, then normalise (the CDS file is wrapped).
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        urllib.request.urlretrieve(NCBI_BASE_URL + NCBI_CDS_FILE, tmp_path)
        with tmp_path.open() as src:
            n = _normalise_fasta(src, fasta)
    finally:
        tmp_path.unlink(missing_ok=True)

    _save_meta("ncbi", {"version": version})
    LOG.info(f"Installed {n} NCBI target genes to {fasta}")
    return fasta


def update_database(database: str, force: bool = False) -> Path:
    """Update the given reference database, returning the installed FASTA path."""
    if database == "card":
        return update_card(force)
    if database == "ncbi":
        return update_ncbi(force)
    raise ValueError(f"Unknown database {database!r}; choose from {DATABASES}")


def _find_member(tf: tarfile.TarFile, suffix: str) -> Optional[tarfile.TarInfo]:
    """Return the first archive member whose name ends with ``suffix`` (or None)."""
    for member in tf.getmembers():
        if member.name.endswith(suffix):
            return member
    return None
