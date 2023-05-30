import os
import sys
from pathlib import Path


def get_pkg_card_fasta_path() -> Path:
    """
    Returns the path to the CARD fasta file in the package.

    https://docs.python.org/3/library/pkgutil.html#pkgutil.get_data
    """
    module_path = Path(os.path.dirname(sys.modules['sarand'].__file__))
    file_path = module_path / 'data' / 'CARD_AMR_seq.fasta'
    return file_path
