import csv
import re
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Optional, List

from sarand.config import PROGRAM_VERSION_NA, CONDA_EXE_NAME, \
    CONDA_RGI_NAME
from sarand.util.logger import LOG


class RgiParams:
    __slots__ = (
        'input_sequence', 'output_file', 'input_type', 'alignment_tool',
        'threads', 'include_loose', 'include_nudge', 'local', 'clean', 'keep',
        'debug', 'low_quality', 'data', 'orf_finder', 'split_prodigal_jobs'
    )

    def __init__(
            self,
            input_sequence: Path,
            output_file: Path,
            input_type: Optional[str] = None,
            alignment_tool: Optional[str] = None,
            threads: Optional[int] = None,
            include_loose: Optional[bool] = None,
            include_nudge: Optional[bool] = None,
            local: Optional[bool] = None,
            clean: Optional[bool] = None,
            keep: Optional[bool] = None,
            debug: Optional[bool] = None,
            low_quality: Optional[bool] = None,
            data: Optional[str] = None,
            orf_finder: Optional[str] = None,
            split_prodigal_jobs: Optional[bool] = None,
    ):
        self.input_sequence: Path = input_sequence
        self.output_file: Path = output_file
        self.input_type: Optional[str] = input_type
        self.alignment_tool: Optional[str] = alignment_tool
        self.threads: Optional[int] = threads
        self.include_loose: Optional[bool] = include_loose
        self.include_nudge: Optional[bool] = include_nudge
        self.local: Optional[bool] = local
        self.clean: Optional[bool] = clean
        self.keep: Optional[bool] = keep
        self.debug: Optional[bool] = debug
        self.low_quality: Optional[bool] = low_quality
        self.data: Optional[str] = data
        self.orf_finder: Optional[str] = orf_finder
        self.split_prodigal_jobs: Optional[bool] = split_prodigal_jobs

    def as_cmd(self) -> List[str]:

        # The version needs to be checked prior to running this, as RGI
        # changed parameters in version 6.
        version = Rgi.version()
        use_v5_params = version[0] in {'4', '5'}

        out = [
            'rgi', 'main',
            '--input_sequence', str(self.input_sequence.absolute()),
            '--output_file', str(self.output_file.absolute()),
        ]
        if self.input_type:
            out.extend(['--input_type', self.input_type])
        if self.alignment_tool:
            out.extend(['--alignment_tool', self.alignment_tool])
        if self.threads:
            out.extend(['--threads', str(self.threads)])
        if self.include_loose:
            out.append('--include_loose')
        if self.include_nudge is not None:
            if use_v5_params:
                if self.include_nudge is False:
                    out.append('--exclude_nudge')
            else:
                if self.include_nudge is True:
                    out.append('--include_nudge')
        if self.local:
            out.append('--local')
        if self.clean:
            out.append('--clean')
        if self.keep:
            out.append('--keep')
        if self.debug:
            out.append('--debug')
        if self.low_quality:
            out.append('--low_quality')
        if self.data:
            out.extend(['--data', self.data])
        if self.orf_finder:
            out.extend(['--orf_finder', self.orf_finder])
        if self.split_prodigal_jobs:
            out.append('--split_prodigal_jobs')
        return out


class RgiResult:
    __slots__ = ('data',)

    def __init__(self, path: Path):
        self.data = self.read_txt(path)

    @staticmethod
    def read_txt(path: Path) -> List[dict]:
        data = list()
        with path.open(newline="") as rgi_file:
            rgi_reader = csv.reader(rgi_file, delimiter="\t")
            next(rgi_reader)
            for row in rgi_reader:
                seq_info = {
                    "ORF_ID": row[0],
                    "gene": row[8].strip(),
                    "prediction_type": row[5].strip(),
                    "best_identities": float(row[9]),
                    "family": row[16].strip(),
                }
                data.append(seq_info)
        return data


class Rgi:
    __slots__ = ('params', 'result')

    def __init__(self, params: RgiParams, result: RgiResult):
        self.params: RgiParams = params
        self.result: RgiResult = result

    @classmethod
    def run(cls, params: RgiParams) -> 'Rgi':
        cmd = params.as_cmd()

        # If this is being run in the Docker container, then activate the env first
        if CONDA_RGI_NAME:
            cmd = [CONDA_EXE_NAME, 'run', '-n', CONDA_RGI_NAME] + cmd

        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()

        # RGI always returns exit code 0...
        if proc.returncode != 0 or len(stderr) > 0:
            LOG.error(stdout)
            LOG.error(stderr)
            raise Exception("ERROR: RGI didn't run successfully!")

        result = RgiResult(Path(f'{params.output_file.absolute()}.txt'))
        return cls(params=params, result=result)

    @classmethod
    def run_for_sarand(cls, input_sequence: Path, output_file: Path, include_loose: bool):
        params = RgiParams(
            input_sequence=input_sequence,
            output_file=output_file,
            input_type='protein',
            clean=True,
            include_nudge=False,
            include_loose=include_loose,
        )
        return cls.run(params=params)

    @staticmethod
    @lru_cache(maxsize=1)
    def version() -> str:
        cmd = ['rgi', '-h']
        if CONDA_RGI_NAME:
            cmd = [CONDA_EXE_NAME, 'run', '-n', CONDA_RGI_NAME] + cmd
        proc = subprocess.Popen(cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            LOG.error(stdout)
            LOG.error(stderr)
            raise Exception("rgi binary not found")
        hits = re.findall(r'(\d\.\d\.\d)', stdout)
        if len(hits) == 1:
            return hits[0]
        return PROGRAM_VERSION_NA
