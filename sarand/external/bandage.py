import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, List, Union, Dict

from sarand.config import PROGRAM_VERSION_NA, CONDA_BANDAGE_NAME, CONDA_EXE_NAME
from sarand.util.logger import LOG

##???????
_RE_VERSION = re.compile(r'(\d+\.\d+\.\d+[-\w]*)')
#????????
_RE_RECORD_CLIP = re.compile(r'(.*)_(\d+)_(\d+)_(\d+)$')


class BandageParams:
    """Parameters for the Bandage software."""

    __slots__ = (
        'graph', 'reads', 'outputfile','pathnodes',
        'minpatcov', 'minhitcov','minmeanid', 'maxevprod',
        'minpatlen', 'maxpatlen','minlendis', 'maxlendis',
        #'threads',
        'verbose'
    )

    def __init__(
            self,
            graph: Optional[Path] = None,
            reads: Optional[Path] = None,
            outputfile: Optional[Path] = None,

            pathnodes: Optional[int] = None,
            minpatcov: Optional[float] = None,
            minmeanid: Optional[float] = None,
            minhitcov: Optional[float] = None,

            #threads: Optional[int] = None,
            verbose: Optional[bool] = None,
            # maxevprod was noted from  sci type in
            # https://manpages.ubuntu.com/manpages/impish/man1/Bandage.1.html
            maxevprod: Optional[float] = None,
            minpatlen: Optional[float] = None,
            maxpatlen: Optional[float] = None,
            minlendis: Optional[int] = None,
            maxlendis: Optional[int] = None,
    ):
        """
                Args:
                    graph: input graph (.gfa)
                    reads: input reads (fasta or fastq, uncompressed or gzipped)
                    outputfile: arg output alignment file (.tsv)

                    pathnodes: The number of allowed nodes in a BLAST query path (1 to 50, default: 6)
                    minpatcov: Minimum fraction of a BLAST query which must be covered by a query path (0.3 to  1, default: 0.9)
                    minmeanid: Minimum mean identity of BLAST hits in a query path (0 to 1, default: 0.5)
                    minhitcov: Minimum  fraction  of  a BLAST query which must be covered by BLAST hits in a query path (0.3 to 1, default: 0.9)

                    threads: number of threads to use
                    verbose: print progress messages
                    maxevprod: Maximum e-value product for all BLAST hits
                        in  a  query  path  (1e-999  to  9.9e1, default: 1e-10)
                    minpatlen: Minimum allowed relative path length as compared
                        to the query (0 to 10000, default:0.95)
                    maxpatlen: Maximum allowed relative path length as compared
                        to the query (0 to 10000, default:1.05)
                    minlendis: Minimum allowed length discrepancy (in bases)
                        between a BLAST query and its path in the graph
                        (-1e+6 to 1e+6, default: off)
                    maxlendis: Maximum allowed length discrepancy (in bases)
                        between a BLAST query and its path in the graph
                        (-1e+6 to 1e+6, default: off)
                """

        # Store the parameters
        self.graph = graph
        self.reads = reads
        self.outputfile = outputfile
        self.pathnodes = pathnodes
        self.minpatcov = minpatcov
        self.minmeanid = minmeanid
        self.minhitcov = minhitcov
        #self.threads = threads
        self.verbose = verbose
        self.maxevprod = maxevprod
        self.minpatlen = minpatlen
        self.maxpatlen = maxpatlen
        self.minlendis = minlendis
        self.maxlendis = maxlendis

    def as_cmd(self) -> List[str]:
        """Return the Bandage command for use in subprocess."""
        cmd = [
            'Bandage',
            'querypaths',
            self.graph.absolute(),
            self.reads.absolute(),
        ]
        if self.outputfile:
            cmd += [self.outputfile.absolute()]
        if self.pathnodes:
            cmd += ['--pathnodes', str(self.pathnodes)]
        if self.minpatcov:
            cmd += ['--minpatcov', str(self.minpatcov)]
        if self.minmeanid:
            cmd += ['--minmeanid', str(self.minmeanid)]
        if self.minhitcov:
            cmd += ['--minhitcov', str(self.minhitcov)]
        #if self.threads:
        #    cmd += ['--threads', str(self.threads)]
        if self.verbose:
            cmd += ['--verbose']
        if self.maxevprod:
            cmd += ['--maxevprod', str(self.maxevprod)]
        if self.minpatlen:
            cmd += ['--minpatlen', str(self.minpatlen)]
        if self.maxpatlen:
            cmd += ['--maxpatlen', str(self.maxpatlen)]
        if self.minlendis:
            cmd += ['--minlendis', str(self.minlendis)]
        if self.maxlendis:
            cmd += ['--maxlendis', str(self.maxlendis)]
        return cmd

    @classmethod
    def from_cli_args(cls, ga):
        """Additional parameters supplied from the command line arguments."""
        d_params = dict()
        if ga is not None:
            for item in ga:
                if len(item) == 2:
                    d_params[item[0]] = item[1]
                elif len(item) == 1:
                    d_params[item[0]] = True
        out = cls()
        out.update_from_dictionary(d_params)
        return out

    def update_from_object(self, other: 'BandageParams'):
        """
        Merge two bandage paramters together, taking the non-None
        values from the other object.
        """
        for k in other.__slots__:
            v = getattr(other, k)
            if v is not None:
                LOG.debug(f'Updated BandageParams.{k} to {v}')
                setattr(self, k, v)
        return

    def update_from_dictionary(self, d):
        """Add/override the parameters from those supplied in a dictionary."""
        if 'graph' in d:
            self.graph = Path(d['graph'])
        if 'reads' in d:
            self.reads = Path(d['reads'])
        if 'outputfile' in d:
            self.outputfile = Path(d['outputfile'])
        if 'pathnodes' in d:
            self.pathnodes = int(d['pathnodes'])
        if 'minpatcov' in d:
            self.minpatcov = float(d['minpatcov'])
        if 'minmeanid' in d:
            self.minmeanid = float(d['minmeanid'])
        if 'minhitcov' in d:
            self.minhitcov = float(d['minhitcov'])
        #if 'threads' in d:
        #    self.threads = int(d['threads'])
        if 'verbose' in d:
            self.verbose = True
        if 'maxevprod' in d:
            self.maxevprod = float(d['maxevprod'])
        if 'minpatlen' in d:
            self.minpatlen = float(d['minpatlen'])
        if 'maxpatlen' in d:
            self.maxpatlen = float(d['maxpatlen'])
        if 'minlendis' in d:
            self.minlendis = int(d['minlendis'])
        if 'maxlendis' in d:
            self.maxlendis = int(d['maxlendis'])
        return


class BandageResult:
    """
    A single entry in the Bandage tsv file.
    """
    __slots__ = (
        'query', 'path_with_start_end', 'length',
        'query_covered_by_path', 'query_covered_by_hits',
        'mean_hit_identity', 'total_hit_mismatches', 'total_hit_gap_opens',
        'relative_length', 'length_discrepancy', 'e_value_product','sequence'
    )

    def __init__(
            self,
            query: str,
            path_with_start_end: str,
            length: int,
            # float but adds %
            query_covered_by_path: str,
            # float but adds %
            query_covered_by_hits: str,
            # float but adds %
            mean_hit_identity: str,
            total_hit_mismatches: int,
            total_hit_gap_opens: int,
            # float but adds %
            relative_length: str,
            length_discrepancy: int,
            e_value_product: float,
            sequence: str
    ):
        self.query: str = query
        self.path_with_start_end: str = path_with_start_end
        self.length: int = length
        # float but adds %
        self.query_covered_by_path: str = query_covered_by_path
        # float but adds %
        self.query_covered_by_hits: str = query_covered_by_hits
        # float but adds %
        self.mean_hit_identity: str = mean_hit_identity
        self.total_hit_mismatches: int = total_hit_mismatches
        self.total_hit_gap_opens: int = total_hit_gap_opens
        # float but adds %
        self.relative_length: str = relative_length
        self.length_discrepancy: int = length_discrepancy
        self.e_value_product: float = e_value_product
        self.sequence: str = sequence

    def __repr__(self):
        return self.query

    @property
    def identity(self) -> str:
        return self.query
        #return self.name.split(' ')[0]
    @property
    def path(self) -> str:
        return (re.sub("\((.*?)\)", "", self.path_with_start_end)).strip()

    @property
    def path_start(self) -> int:
        if self.path_with_start_end.startswith("("):
            index = self.path_with_start_end.find(")")
            return int(self.path_with_start_end[1 : index])
        return 0

    @property
    def path_end(self) -> int:
        if self.path_with_start_end.endswith(")"):
            index = self.path_with_start_end.rfind("(")
            return int(self.path_with_start_end[index + 1 : -1])
        return 0

    @property
    def path_to_sarand(self):
        """
        self.path: tsv files from Bandage displays a list of nodes
            with -/+ tail and comma separated : e.g., '(1363) 69-, 2193+ (1786)'
        Return:
            nodes:	list of node numbers -> e.g., [69, 2193]
            orientation: list of orientation of nodes -> e.g., [-, +]
        """
        # Remove text between ()
        #purePath = (re.sub("\((.*?)\)", "", self.path_with_start_end)).strip()
        if len(self.path) == 0:
            raise Exception('??')
        nodes = list()
        orientation = list()
        node_list = self.path.split(",")
        for node in node_list:
            if "-" in node:
                orientation.append("-")
            else:
                orientation.append("+")
            node = re.sub("[+-]", "", node.split()[0])
            nodes.append(node)
        return nodes, orientation

    @property
    def coverage_pct(self) -> float:
        """Calculate the percentage coverage of the alignment."""
        return float(re.sub("[%]", "", self.query_covered_by_path))

    @property
    def identity_pct(self) -> float:
        """Calculate the percentage identity of the alignment."""
        return float(re.sub("[%]", "", self.mean_hit_identity))

    @property
    def amr_name(self) -> str:
        amr_str = self.identity.split("|")[-1].strip().replace(" ", "_").replace("'", ";")
        return BandageResult.restricted_amr_name_from_modified_name(amr_str)

    @staticmethod
    def restricted_amr_name_from_modified_name(amr_name):
        """Lifted from utils to avoid circular imports
        TODO: Refactor utils to avoid calls to Bandage?
        """
        amr_name1 = amr_name.replace(";", "SS")
        amr_name1 = "".join(
            e for e in amr_name1 if e.isalpha() or e.isnumeric() or e == "_" or e == "-"
        )
        return amr_name1


class Bandage:
    """Wrapper to the Bandage software.
    https://github.com/rrwick/Bandage
    """
    __slots__ = ('params', 'results', 'stdout', 'stderr')

    def __init__(
            self,
            params: BandageParams,
            results: List[BandageResult],
            stdout: Optional[str] = None,
            stderr: Optional[str] = None,
    ):
        self.params: BandageParams = params
        self.results: List[BandageResult] = results
        self.stdout: Optional[str] = stdout
        self.stderr: Optional[str] = stderr

    @classmethod
    def run(cls, params: BandageParams) -> 'Bandage':
        """Runs Bandage with the given parameters."""

        # Display the command to be run
        cmd = params.as_cmd()

        # If this is being run in the Docker container, then activate the env first
        if CONDA_BANDAGE_NAME:
            cmd = [CONDA_EXE_NAME, 'run', '-n', CONDA_BANDAGE_NAME] + cmd
        LOG.info(' '.join(map(str, cmd)))

        # Run the command
        proc = subprocess.Popen(cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            LOG.error(stdout)
            LOG.error(stderr)
            raise Exception('Bandage binary not found')

        results = Bandage.read_file(params.outputfile.with_suffix('.tsv'))
        return cls(params, results, stdout, stderr)

    @classmethod
    def run_for_sarand(
            cls,
            gfa: Path,
            reads: Path,
            threshold: float,
            ga_extra_args: BandageParams,
            out_dir: Optional[Path] = None,
            #threads: Optional[int] = 1,
    ) -> 'Bandage':

        """Default method to run Bandage for the sarand pipeline."""
        if isinstance(reads, str):
            reads = Path(reads)

        if out_dir is not None:
            params = BandageParams(
                graph=gfa,
                reads=reads,
                outputfile=out_dir / 'bandage',
                minpatcov = ((threshold - 1) / 100.0),
                minmeanid = ((threshold - 1) / 100.0),
                minhitcov = ((threshold - 1) / 100.0),
                #threads=threads,
            )
            if ga_extra_args:
                params.update_from_object(ga_extra_args)
            out = Bandage.run(params)
            if out.stdout:
                path_stdout = out_dir / 'bandage_stdout.txt'
                with path_stdout.open('w') as f:
                    f.write(out.stdout)
            if out.stderr:
                path_stderr = out_dir / 'bandage_stderr.txt'
                with path_stderr.open('w') as f:
                    f.write(out.stderr)
            return out

        # Run Bandage in a temporary directory
        else:
            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_dir = Path(tmp_dir)
                params = BandageParams(
                    graph=gfa,
                    reads=reads,
                    outputfile=tmp_dir / 'bandage',
                    minpatcov = ((threshold - 1) / 100.0),
                    minmeanid = ((threshold - 1) / 100.0),
                    minhitcov = ((threshold - 1) / 100.0),
                    #threads=threads,
                )
                if ga_extra_args:
                    params.update_from_object(ga_extra_args)
                return Bandage.run(params)

    @staticmethod
    def read_file(path: Path) -> List[BandageResult]:

        out = list()
        with path.open() as f:
            for line in f.readlines()[1:]:
                cols = line.strip().split('\t')
                query = cols[0]
                path_with_start_end = cols[1]
                length = int(cols[2])
                query_covered_by_path = cols[3]
                query_covered_by_hits = cols[4]
                mean_hit_identity = cols[5]
                total_hit_mismatches = int(cols[6])
                total_hit_gap_opens = int(cols[7])
                relative_length = cols[8]
                length_discrepancy = int(cols[9])
                e_value_product = float(cols[10])
                sequence = cols[11]

                obj = BandageResult(
                    query, path_with_start_end, length, query_covered_by_path,
                    query_covered_by_hits, mean_hit_identity,
                    total_hit_mismatches, total_hit_gap_opens, relative_length,
                    length_discrepancy, e_value_product, sequence
                )
                out.append(obj)
        return out

    @staticmethod
    def version() -> str:
        """Returns the version of Bandage on the path."""
        cmd = ['Bandage', '--version']
        if CONDA_BANDAGE_NAME:
            cmd = [CONDA_EXE_NAME, 'run', '-n', CONDA_BANDAGE_NAME] + cmd
        LOG.debug(' '.join(map(str, cmd)))
        proc = subprocess.Popen(
            cmd,
            encoding='utf-8',
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise Exception('Bandage binary not found')
        match_stdout, match_stderr = None, None
        if stdout:
            match_stdout = _RE_VERSION.search(stdout)
        if stderr:
            match_stderr = _RE_VERSION.search(stderr)
        for match in (match_stdout, match_stderr):
            if match:
                ver = str(match.group(1))
                if ver.endswith('-'):
                    ver = ver[:-1]
                return ver
        return PROGRAM_VERSION_NA
