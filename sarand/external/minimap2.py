"""Wrapper for the minimap2 aligner (used for batched sequence comparison).

sarand uses minimap2 to decide whether ORF nucleotide sequences are
near-identical when collapsing duplicate neighborhoods. ORFs are short, so the
seeding is tuned (``-k 11 -w 5``) to map sequences down to ~90 bp, which the
default/asm presets miss; ``-c`` requests base-level alignment so the reported
match/block lengths give an accurate percent identity.
"""
import re
import subprocess
from pathlib import Path
from typing import List

from sarand.config import PROGRAM_VERSION_NA
from sarand.util.logger import LOG

# minimap2 --version prints e.g. "2.31-r1302"
_RE_VERSION = re.compile(r'(\d+\.\d+[.\-\w]*)')

# Seeding tuned for short, high-identity ORF nucleotide sequences.
_SHORT_SEQ_ARGS = ['-c', '-k', '11', '-w', '5']


class Minimap2Result:
    """A single PAF record from a minimap2 run."""

    __slots__ = (
        'query', 'query_len', 'target', 'target_len', 'num_matches', 'block_len'
    )

    def __init__(self, query: str, query_len: int, target: str, target_len: int,
                 num_matches: int, block_len: int):
        """Store the PAF columns needed to derive identity and coverage."""
        self.query: str = query
        self.query_len: int = query_len
        self.target: str = target
        self.target_len: int = target_len
        self.num_matches: int = num_matches
        self.block_len: int = block_len

    @property
    def identity_pct(self) -> float:
        """Percent identity over the aligned block (matches / block length)."""
        if self.block_len == 0:
            return 0.0
        return self.num_matches / self.block_len * 100

    @property
    def coverage_pct(self) -> float:
        """Percent of the longer sequence covered by the aligned block."""
        longer = max(self.query_len, self.target_len)
        if longer == 0:
            return 0.0
        return self.block_len / longer * 100


class Minimap2:
    """Wrapper for the minimap2 program."""

    @staticmethod
    def run(query: Path, target: Path) -> List['Minimap2Result']:
        """Align all sequences in ``query`` against all in ``target`` (one call).

        Parameters:
            query: FASTA of query sequences
            target: FASTA of reference sequences (indexed by minimap2)
        Return:
            the parsed PAF records (one or more per aligned query/target pair)
        """
        cmd = [
            'minimap2', *_SHORT_SEQ_ARGS,
            str(Path(target).absolute()),
            str(Path(query).absolute()),
        ]
        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8'
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f'minimap2 failed: {stderr}')
        return Minimap2.parse_paf(stdout)

    @staticmethod
    def parse_paf(data: str) -> List['Minimap2Result']:
        """Parse minimap2 PAF stdout into a list of results."""
        out = []
        for line in data.splitlines():
            cols = line.split('\t')
            if len(cols) < 11:
                continue
            out.append(Minimap2Result(
                query=cols[0],
                query_len=int(cols[1]),
                target=cols[5],
                target_len=int(cols[6]),
                num_matches=int(cols[9]),
                block_len=int(cols[10]),
            ))
        return out

    @staticmethod
    def version() -> str:
        """Returns the version of minimap2 on the path."""
        cmd = ['minimap2', '--version']
        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(
            cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise Exception('minimap2 binary not found')
        re_hit = _RE_VERSION.search(stdout)
        if re_hit:
            return re_hit.group(1)
        return PROGRAM_VERSION_NA
