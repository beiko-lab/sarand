import re
import subprocess
from enum import Enum
from pathlib import Path
from typing import List, Optional

from sarand.config import PROGRAM_VERSION_NA
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
            identity: int,
            word_length: int,
            threads: int,
    ):
        self.input_file: Path = input_file
        self.output_file: Path = output_file 
        self.identitiy: int = identity 
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
            self.identity,
            '-n',
            self.word_length,
            '-T',
            self.threads
        ]
        return out


class Cdhit:
    """Wrapper for the cd-hit program."""

    __slots__ = ('params')

    def __init__(self, params: CdhitParams):
        self.params: CdhitParams = params

    @classmethod
    def run(cls, params: BlastnParams):
        """Run cd_hit given the parameters."""

        # Generate the command
        cmd = params.as_cmd()

        # Run blastn
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

    @staticmethod
    def version() -> str:
        """Returns the version of blastn on the path."""
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
        else:
            raise Exception('valid cd-hit binary not found')
        return PROGRAM_VERSION_NA
