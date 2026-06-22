"""Unit tests for sarand.util.cli argparse type helpers."""
import argparse

import pytest

from sarand.util.cli import check_file, validate_range


def test_check_file_returns_resolved_path(tmp_path):
    f = tmp_path / "input.gfa"
    f.write_text("data")
    result = check_file(str(f))
    assert result == f.resolve()
    assert result.is_absolute()


def test_check_file_missing_raises(tmp_path):
    with pytest.raises(argparse.ArgumentTypeError):
        check_file(str(tmp_path / "does_not_exist.gfa"))


def test_check_file_directory_raises(tmp_path):
    """A directory is not a valid input file."""
    with pytest.raises(argparse.ArgumentTypeError):
        check_file(str(tmp_path))


def test_validate_range_float_in_range():
    checker = validate_range(float, 0.0, 1.0)
    assert checker("0.5") == 0.5
    assert checker("0") == 0.0
    assert checker("1") == 1.0


def test_validate_range_int_in_range():
    checker = validate_range(int, 1, 50)
    assert checker("6") == 6
    assert isinstance(checker("6"), int)


@pytest.mark.parametrize("arg", ["-0.1", "1.1", "100"])
def test_validate_range_out_of_range_raises(arg):
    checker = validate_range(float, 0.0, 1.0)
    with pytest.raises(argparse.ArgumentTypeError):
        checker(arg)


def test_validate_range_float_unparseable_raises():
    checker = validate_range(float, 0.0, 1.0)
    with pytest.raises(argparse.ArgumentTypeError):
        checker("not_a_number")


def test_validate_range_int_unparseable_raises():
    checker = validate_range(int, 1, 50)
    with pytest.raises(argparse.ArgumentTypeError):
        checker("3.5")
