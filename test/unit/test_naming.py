"""Unit tests for sarand.util.naming string helpers."""
from pathlib import Path

import pytest

from sarand.util.naming import (
    extract_name_from_file_name,
    restricted_target_name_from_modified_name,
    target_name_from_comment,
)


@pytest.mark.parametrize(
    "file_name,expected",
    [
        ("ErmA.fasta", "ErmA"),
        ("/some/dir/ANT9-Ia.txt", "ANT9-Ia"),
        (Path("/a/b/c.fasta"), "c"),
        ("no_extension", "no_extension"),
        ("archive.tar.gz", "archive.tar"),
    ],
)
def test_extract_name_from_file_name(file_name, expected):
    assert extract_name_from_file_name(file_name) == expected


@pytest.mark.parametrize(
    "comment,expected",
    [
        # the gene id is the last pipe-delimited field of the first
        # whitespace-delimited token (the "[Staphylococcus]" part is dropped)
        ("gb|X|ARO|aac(6')-Ie [Staphylococcus]", "aac(6;)-Ie"),
        # only the header up to the first whitespace is used (matching how
        # BLAST/Bandage truncate query names), so a trailing token is dropped
        ("gb|X|ARO|some gene name", "some"),
        # NCBI AMR_CDS.fa headers carry a trailing "accession:coords" token
        # after the gene name; it must be dropped so the name matches Bandage
        ("AAA16360.1|L11078.1|1|1|stxA2b|stxA2b|stxA2b L11078.1:177-1136", "stxA2b"),
        # forward slash becomes ']'
        ("foo/bar", "foo]bar"),
        # apostrophe becomes ';'
        ("aac(6')", "aac(6;)"),
        # surrounding whitespace is stripped
        ("   trimmed   ", "trimmed"),
    ],
)
def test_target_name_from_comment(comment, expected):
    assert target_name_from_comment(comment) == expected


@pytest.mark.parametrize(
    "name,expected",
    [
        # ';' expands to 'SS', parentheses are dropped
        ("aac(6;)-Ie", "aac6SS-Ie"),
        # underscores and hyphens survive, other punctuation removed
        ("ANT9-Ia", "ANT9-Ia"),
        ("a_b-c", "a_b-c"),
        ("name with spaces!", "namewithspaces"),
        ("keep123", "keep123"),
    ],
)
def test_restricted_target_name_from_modified_name(name, expected):
    assert restricted_target_name_from_modified_name(name) == expected
