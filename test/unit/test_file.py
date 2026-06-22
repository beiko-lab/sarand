"""Unit tests for sarand.util.file."""
import json

from sarand.util.file import try_dump_to_disk


def test_try_dump_to_disk_writes_json(tmp_path):
    obj = {"a": 1, "b": [1, 2, 3]}
    path = tmp_path / "out.json"
    try_dump_to_disk(obj, path)
    assert json.loads(path.read_text()) == obj


def test_try_dump_to_disk_swallows_errors(tmp_path):
    """A non-serialisable object must not raise (best-effort debug dump)."""
    path = tmp_path / "out.json"
    # sets are not JSON serialisable; the helper should log and return cleanly
    try_dump_to_disk({1, 2, 3}, path)
