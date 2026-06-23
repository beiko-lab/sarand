import re
import subprocess
from pathlib import Path
from typing import List

from sarand.util.logger import LOG

_RE_VERSION = re.compile(r'		====== CD-HIT version ([0-9]+\.[0-9]+\.[0-9]+) \(built on .*\) ======')

class CdhitParams:
    """Wrapper for cd-hit parameters that have been implemented."""

    __slots__ = (
        'input_file', 'output_file', 'identity', 'word_length', 'threads'
    )

    def __init__(
            self,
            input_file: Path,
            output_file: Path,
            identity: float,
            word_length: int,
            threads: int,
    ):
        """Store the cd-hit parameters (input/output paths, identity, word length, threads)."""
        self.input_file: Path = input_file
        self.output_file: Path = output_file
        self.identity: float = identity
        self.word_length: int = word_length
        self.threads: int = threads

    def as_cmd(self) -> List[str]:
        """Generate the command as a list of strings."""

        out = [
            'cd-hit',
            '-i',
            str(self.input_file.absolute()),
            '-o',
            str(self.output_file.absolute()),
            '-c',
            str(self.identity),
            '-n',
            str(self.word_length),
            '-T',
            str(self.threads),
            '-M', # unlimited memory for cd-hit instead of built in default limit
            '0'
        ]
        return out


class Cdhit:
    """Wrapper for the cd-hit program."""

    __slots__ = ('params')

    def __init__(self, params: CdhitParams):
        """Wrap a completed cd-hit run described by ``params``."""
        self.params: CdhitParams = params

    @classmethod
    def run(cls, params: CdhitParams):
        """Run cd-hit given the parameters."""

        # Generate the command
        cmd = params.as_cmd()

        # Run cd-hit
        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding='utf-8'
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f'cd-hit failed: {stderr}')

        return cls(params)

    @classmethod
    def cluster(
            cls,
            input_file,
            output_file,
            identity: float,
            word_length: int = 5,
            threads: int = 1,
    ):
        """Cluster the sequences in input_file at the given identity.

        word_length defaults to cd-hit's default of 5 and threads to 1, since
        per-target parallelism is handled by the calling pool.
        """
        params = CdhitParams(
            input_file=Path(input_file),
            output_file=Path(output_file),
            identity=identity,
            word_length=word_length,
            threads=threads,
        )
        return cls.run(params)

    @staticmethod
    def version() -> str:
        """Return the version of the cd-hit binary on the path."""
        cmd = ['cd-hit', '-h']
        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(
            cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        # cd-hit always returns 1 unless actually executed so can't use this flag
        re_hit = _RE_VERSION.match(stdout)
        if re_hit:
            return re_hit.group(1)
        raise Exception('valid cd-hit binary not found')
