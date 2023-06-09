import json
from pathlib import Path

from sarand.util.logger import LOG


def try_dump_to_disk(obj, path: Path):
    """Dump a JSON serializable object to disk (for debugging)"""
    LOG.debug(f'Dumping object to disk: {path.absolute()}')
    try:
        with open(path, 'w') as f:
            json.dump(obj, f, indent=2)
    except Exception as e:
        LOG.error(f'Failed to dump object to disk: {e}')
    return
