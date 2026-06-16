import json
from pathlib import Path
from typing import Any

from sarand.util.logger import LOG


def try_dump_to_disk(obj: Any, path: Path) -> None:
    """Dump a JSON-serializable object to ``path`` (best-effort; used for debugging)."""
    LOG.debug(f'Dumping object to disk: {path.absolute()}')
    try:
        with open(path, 'w') as f:
            json.dump(obj, f, indent=2)
    except Exception as e:
        LOG.error(f'Failed to dump object to disk: {e}')
    return
