import re
import subprocess
from enum import Enum
from pathlib import Path
from typing import List, Optional

from sarand.config import PROGRAM_VERSION_NA, CONDA_BLAST_NAME, CONDA_EXE_NAME
from sarand.util.logger import LOG

_RE_VERSION = re.compile(r'blastn: ([\d.]+)')


class BlastnOutFmt(Enum):
    """Enumerable class of all supported output formats."""
    FMT_1 = "10 pident length qcovhsp"


class BlastnParams:
    """Wrapper for blastn parameters that have been implemented."""

    __slots__ = (
        'query', 'subject', 'task', 'outfmt', 'max_target_seqs', 'evalue',
        'perc_identity'
    )

    def __init__(
            self,
            query: Path,
            subject: Path,
            task: str,
            outfmt: BlastnOutFmt,
            max_target_seqs: Optional[int] = None,
            evalue: Optional[float] = None,
            perc_identity: Optional[float] = None,
    ):
        self.query: Path = query
        self.subject: Path = subject
        self.task: str = task
        self.outfmt: BlastnOutFmt = outfmt
        self.max_target_seqs: Optional[int] = max_target_seqs
        self.evalue: Optional[float] = evalue
        self.perc_identity: Optional[float] = perc_identity

    def as_cmd(self) -> List[str]:
        """Generate the command as a list of strings."""

        out = [
            'blastn',
            '-query',
            str(self.query.absolute()),
            '-subject',
            str(self.subject.absolute()),
            '-task',
            self.task,
            '-outfmt',
            str(self.outfmt.value),
        ]
        if self.evalue is not None:
            out += ['-evalue', str(self.evalue)]
        if self.perc_identity is not None:
            out += ['-perc_identity', str(self.perc_identity)]
        if self.max_target_seqs is not None:
            out += ['-max_target_seqs', str(self.max_target_seqs)]
        return out


class BlastnResult:
    """Wrapper for a single blastn result."""

    __slots__ = ('pident', 'length', 'qcovhsp')

    def __init__(
            self,
            pident: float,
            length: int,
            qcovhsp: int,
    ):
        self.pident: float = pident
        self.length: int = length
        self.qcovhsp: int = qcovhsp

    @classmethod
    def from_outfmt(cls, data: str, outfmt: BlastnOutFmt) -> 'BlastnResult':
        """Parse a single blastn result from the given output format."""

        if outfmt is BlastnOutFmt.FMT_1:
            cols = data.split(',')
            return cls(
                pident=float(cols[0]),
                length=int(cols[1]),
                qcovhsp=int(cols[2]),
            )

        raise NotImplementedError(f'outfmt {outfmt} not implemented')


class Blastn:
    """Wrapper for the blastn program."""

    __slots__ = ('params', 'results')

    def __init__(self, params: BlastnParams, results: List[BlastnResult]):
        self.params: BlastnParams = params
        self.results: List[BlastnResult] = results

    @classmethod
    def run(cls, params: BlastnParams):
        """Run blastn given the parameters."""

        # Generate the command
        cmd = params.as_cmd()

        # Run blastn
        if CONDA_BLAST_NAME:
            cmd = [
                      CONDA_EXE_NAME,
                      'run',
                      '-n',
                      CONDA_BLAST_NAME,
                  ] + cmd
        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding='utf-8'
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f'blastn failed: {stderr}')

        # Parse the output
        results = Blastn.from_outfmt(stdout, params.outfmt)

        # Return the class
        return cls(params, results)

    @classmethod
    def run_for_sarand_compare_two_sequences(cls, query: Path, subject: Path):
        """Special implementation for sarand compare two sequences."""

        params = BlastnParams(
            query=query,
            subject=subject,
            task='blastn-short',
            outfmt=BlastnOutFmt.FMT_1,
        )
        return Blastn.run(params)

    @staticmethod
    def from_outfmt(data: str, outfmt: BlastnOutFmt) -> List[BlastnResult]:
        """Read the stdout of a blastn call and convert it to results."""

        out = list()
        for line in data.splitlines():
            out.append(BlastnResult.from_outfmt(line, outfmt))
        return out

    @staticmethod
    def version() -> str:
        """Returns the version of blastn on the path."""
        cmd = ['blastn', '-version']
        if CONDA_BLAST_NAME:
            cmd = [
                      CONDA_EXE_NAME,
                      'run',
                      '-n',
                      CONDA_BLAST_NAME,
                  ] + cmd
        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(
            cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise Exception('blastn binary not found')
        re_hit = _RE_VERSION.match(stdout)
        if re_hit:
            return re_hit.group(1)
        return PROGRAM_VERSION_NA
