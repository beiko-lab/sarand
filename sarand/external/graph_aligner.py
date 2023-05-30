import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional, List, Union, Dict, Tuple

from Bio import SeqIO

_RE_VERSION = re.compile(r'(\d+\.\d+\.\d+[-\w]*)')
_RE_RECORD_CLIP = re.compile(r'.*_(\d+)_(\d+)_(\d+)$')


class GraphAlignerParams:
    """Parameters for the GraphAligner software."""

    __slots__ = ('graph', 'reads', 'align_out', 'corrected_out', 'corrected_clipped_out',
                 'threads', 'verbose', 'e_cutoff', 'min_align_score', 'multimap_score_fraction',
                 'max_cluster_extend', 'seeds_cluster_size', 'seeds_minimizer_length',
                 'seeds_minimizer_windowsize', 'seeds_minimizer_density',
                 'seeds_minimizer_ignore_frequent', 'seeds_mum_count', 'seeds_mem_count',
                 'seeds_mxm_length', 'seeds_mxm_cache_prefix', 'seeds_mxm_windowsize',
                 'bandwidth', 'tangle_effort', 'x_drop', 'precise_clipping', 'max_trace_count',
                 'preset'
                 )

    def __init__(
            self,
            graph: Path,
            reads: Path,
            align_out: Optional[Path] = None,
            corrected_out: Optional[Path] = None,
            corrected_clipped_out: Optional[Path] = None,

            threads: Optional[int] = None,
            verbose: bool = False,
            e_cutoff: Optional[float] = None,
            min_align_score: Optional[float] = None,
            multimap_score_fraction: Optional[float] = None,

            max_cluster_extend: Optional[int] = None,
            seeds_cluster_size: Optional[int] = None,
            seeds_minimizer_length: Optional[int] = None,
            seeds_minimizer_windowsize: Optional[int] = None,
            seeds_minimizer_density: Optional[float] = None,
            seeds_minimizer_ignore_frequent: Optional[float] = None,
            seeds_mum_count: Optional[int] = None,
            seeds_mem_count: Optional[int] = None,
            seeds_mxm_length: Optional[int] = None,
            seeds_mxm_cache_prefix: Optional[str] = None,
            seeds_mxm_windowsize: Optional[int] = None,

            bandwidth: Optional[int] = None,
            tangle_effort: Optional[int] = None,
            x_drop: Optional[int] = None,
            precise_clipping: Optional[float] = None,
            max_trace_count: Optional[int] = None,

            preset: Optional[str] = None,
    ):
        """
                Args:
                    graph: input graph (.gfa / .vg)
                    reads: input reads (fasta or fastq, uncompressed or gzipped)
                    align_out: arg output alignment file (.gaf/.gam/.json)
                    corrected_out: output corrected reads file (.fa/.fa.gz)
                    corrected_clipped_out: output corrected clipped reads file (.fa/.fa.gz)

                    threads: number of threads to use
                    verbose: print progress messages
                    e_cutoff: discard alignments with E-value > arg
                    min_align_score: discard alignments with alignment score < arg (default 0)
                    discard alignments whose alignment score is less than this fraction of the best overlapping alignment (default 0.9)

                    max_cluster_extend:  extend up to arg seed clusters (-1 for all) (default 10)
                    seeds_cluster_size: discard seed clusters with fewer than arg seeds
                    seeds_minimizer_length: k-mer length for minimizer seeding
                    seeds_minimizer_windowsize: window size for minimizer seeding
                    seeds_minimizer_density: keep approximately (arg * sequence length) least frequent minimizers (-1 for all)
                    seeds_minimizer_ignore_frequent: ignore arg most frequent fraction of minimizers
                    seeds_mum_count: arg longest maximal unique matches (-1 for all)
                    seeds_mem_count: arg longest maximal exact matches (-1 for all)
                    seeds_mxm_length: minimum length for maximal unique / exact matches
                    seeds_mxm_cache_prefix: store the mum/mem seeding index to the disk for reuse, or reuse it if it exists
                    seeds_mxm_windowsize: window size for mem/mum seeding (0 for no windowing)

                    bandwidth: alignment bandwidth
                    tangle_effort: tangle effort limit (-1 for unlimited)
                    x_drop: X-drop alignment ending score cutoff
                    precise_clipping: clip the alignment ends with arg as the identity cutoff between correct / wrong alignments (default 0.66)
                    max_trace_count: backtrace from up to arg highest scoring local maxima per cluster (-1 for all)

                    preset: Parameters optimized for de Bruijn or variation graphs: dbg, or vg respectively
                """
        self.graph = graph
        self.reads = reads
        self.align_out = align_out
        self.corrected_out = corrected_out
        self.corrected_clipped_out = corrected_clipped_out
        self.threads = threads
        self.verbose = verbose
        self.e_cutoff = e_cutoff
        self.min_align_score = min_align_score
        self.multimap_score_fraction = multimap_score_fraction
        self.max_cluster_extend = max_cluster_extend
        self.seeds_cluster_size = seeds_cluster_size
        self.seeds_minimizer_length = seeds_minimizer_length
        self.seeds_minimizer_windowsize = seeds_minimizer_windowsize
        self.seeds_minimizer_density = seeds_minimizer_density
        self.seeds_minimizer_ignore_frequent = seeds_minimizer_ignore_frequent
        self.seeds_mum_count = seeds_mum_count
        self.seeds_mem_count = seeds_mem_count
        self.seeds_mxm_length = seeds_mxm_length
        self.seeds_mxm_cache_prefix = seeds_mxm_cache_prefix
        self.seeds_mxm_windowsize = seeds_mxm_windowsize
        self.bandwidth = bandwidth
        self.tangle_effort = tangle_effort
        self.x_drop = x_drop
        self.precise_clipping = precise_clipping
        self.max_trace_count = max_trace_count
        self.preset = preset

    def as_cmd(self) -> List[str]:
        cmd = [
            'GraphAligner',
            '--graph', self.graph.absolute(),
            '--reads', self.reads.absolute(),
        ]
        if self.align_out:
            cmd += ['--alignments-out', self.align_out.absolute()]
        if self.corrected_out:
            cmd += ['--corrected-out', self.corrected_out.absolute()]
        if self.corrected_clipped_out:
            cmd += ['--corrected-clipped-out', self.corrected_clipped_out.absolute()]
        if self.threads:
            cmd += ['--threads', str(self.threads)]
        if self.verbose:
            cmd += ['--verbose']
        if self.e_cutoff:
            cmd += ['--E-cutoff', str(self.e_cutoff)]
        if self.min_align_score:
            cmd += ['--min-alignment-score', str(self.min_align_score)]
        if self.multimap_score_fraction:
            cmd += ['--multimap-score-fraction', str(self.multimap_score_fraction)]
        if self.max_cluster_extend:
            cmd += ['--max-cluster-extend', str(self.max_cluster_extend)]
        if self.seeds_cluster_size:
            cmd += ['--seeds-clustersize', str(self.seeds_cluster_size)]
        if self.seeds_minimizer_length:
            cmd += ['--seeds-minimizer-length', str(self.seeds_minimizer_length)]
        if self.seeds_minimizer_windowsize:
            cmd += ['--seeds-minimizer-windowsize', str(self.seeds_minimizer_windowsize)]
        if self.seeds_minimizer_density:
            cmd += ['--seeds-minimizer-density', str(self.seeds_minimizer_density)]
        if self.seeds_minimizer_ignore_frequent:
            cmd += ['--seeds-minimizer-ignore-frequent', str(self.seeds_minimizer_ignore_frequent)]
        if self.seeds_mum_count:
            cmd += ['--seeds-mum-count', str(self.seeds_mum_count)]
        if self.seeds_mem_count:
            cmd += ['--seeds-mem-count', str(self.seeds_mem_count)]
        if self.seeds_mxm_length:
            cmd += ['--seeds-mxm-length', str(self.seeds_mxm_length)]
        if self.seeds_mxm_cache_prefix:
            cmd += ['--seeds-mxm-cache-prefix', str(self.seeds_mxm_cache_prefix)]
        if self.seeds_mxm_windowsize:
            cmd += ['--seeds-mxm-windowsize', str(self.seeds_mxm_windowsize)]
        if self.bandwidth:
            cmd += ['--bandwidth', str(self.bandwidth)]
        if self.tangle_effort:
            cmd += ['--tangle-effort', str(self.tangle_effort)]
        if self.x_drop:
            cmd += ['--X-drop', str(self.x_drop)]
        if self.precise_clipping:
            cmd += ['--precise-clipping', str(self.precise_clipping)]
        if self.max_trace_count:
            cmd += ['--max-trace-count', str(self.max_trace_count)]
        if self.preset:
            cmd += ['--preset', self.preset]
        return cmd


class GraphAlignerResult:
    """
    https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
    """
    __slots__ = (
        'name', 'len', 'start', 'end', 'strand', 'path', 'path_len',
        'path_start', 'path_end', 'n_res_matches', 'aln_block_len',
        'map_quality', 'custom_tags', 'sequence'
    )

    def __init__(
            self,
            seq_name: str,
            seq_len: int,
            seq_start: int,
            seq_end: int,
            seq_strand: str,
            seq_path: str,
            path_length: int,
            path_start: int,
            path_end: int,
            n_res_matches: int,
            aln_block_len: int,
            map_quality: float,
            custom_tags: Dict[str, Union[float, str, int]],
            sequence: Optional[str],
    ):
        self.name: str = seq_name
        self.len: int = seq_len
        self.start: int = seq_start
        self.end: int = seq_end
        self.strand: str = seq_strand
        self.path: str = seq_path
        self.path_len: int = path_length
        self.path_start: int = path_start
        self.path_end: int = path_end
        self.n_res_matches: int = n_res_matches
        self.aln_block_len: int = aln_block_len
        self.map_quality: float = map_quality
        self.custom_tags: Dict[str, Union[float, str, int]] = custom_tags
        self.sequence: Optional[str] = sequence

        """
        Custom tags:
        https://github.com/maickrau/GraphAligner/blob/master/src/GraphAlignerGAFAlignment.h#LL193C3-L193C3
      sstr << "\t" << "NM:i:" << (mismatches + deletions + insertions);
		if (alignmentXScore != -1) sstr << "\t" << "AS:f:" << alignmentXScore;
		sstr << "\t" << "dv:f:" << 1.0-((double)matches / (double)(matches + mismatches + deletions + insertions));
		sstr << "\t" << "id:f:" << ((double)matches / (double)(matches + mismatches + deletions + insertions));
		if (includecigar) sstr << "\t" << "cg:Z:" << cigar.str();
		return sstr.str();
		
        https://github.com/pangenome/rs-peanut
        """

    def __repr__(self):
        return self.name

    @property
    def identity(self) -> str:
        return self.name.split(' ')[0]

    @property
    def path_to_sarand(self):
        hits = re.findall(r'([><])(\d+)', self.path)
        if len(hits) == 0:
            raise Exception('??')
        nodes = list()
        orientation = list()
        for (strand, node_id) in hits:
            nodes.append(node_id)
            orientation.append('-' if strand == '<' else '+')
        return nodes, orientation

    @property
    def coverage_pct(self) -> float:
        """Calculate the coverage percentage of the alignment.
        TODO: Validate that this is the correct calculation
        """
        return self.n_res_matches / self.len * 100

    @property
    def identity_pct(self) -> float:
        """Calculate the percentage identity of the alignment.
        TODO: Validate that this is the correct calculation
        """
        return self.custom_tags['id'] * 100


class GraphAligner:
    """Wrapper to the GraphAligner software.
    https://github.com/maickrau/GraphAligner
    """

    def __init__(self, params: GraphAlignerParams, results: List[GraphAlignerResult]):
        self.params = params
        self.results: List[GraphAlignerResult] = results

    @classmethod
    def run(cls, params: GraphAlignerParams) -> 'GraphAligner':
        """Runs GraphAligner with the given parameters."""

        # Limit the available options as not all methods have been implemented
        if params.corrected_clipped_out is None or params.corrected_clipped_out.suffix != '.fa':
            raise NotImplemented("corrected_clipped_out must be specified with a .fa suffix")
        if params.align_out is None or params.align_out.suffix != '.gaf':
            raise NotImplemented("align_out must be specified with a .gaf suffix")

        # Display the command to be run
        cmd = params.as_cmd()
        # print(' '.join(map(str, cmd)))

        # Run the command
        proc = subprocess.Popen(cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            print(stdout, stderr)
            raise Exception('GraphAligner binary not found')

        # Process the output
        results = GraphAligner.read_file(params.align_out)
        results_clip = GraphAligner.read_clip_fasta(params.corrected_clipped_out)

        # Add the sequence to each of the results
        for result in results:
            result.sequence = results_clip[result.identity][(result.start, result.end)]
        return cls(params, results)

    @classmethod
    def run_for_sarand(cls, gfa: Path, reads: Path, threshold: float) -> 'GraphAligner':
        """Default method to run GraphAligner for the sarand pipeline."""
        if isinstance(reads, str):
            # print('Converitng reads to path object')
            reads = Path(reads)
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            params = GraphAlignerParams(
                graph=gfa,
                reads=reads,
                align_out=tmp_dir / 'align.gaf',
                corrected_clipped_out=tmp_dir / 'corrected_clipped.fa',
                min_align_score=threshold,
                preset='vg',
            )
            return GraphAligner.run(params)

    @staticmethod
    def read_clip_fasta(path: Path) -> Dict[str, Dict[Tuple[int, int], str]]:
        out = defaultdict(dict)
        with path.open() as f:
            for record in SeqIO.parse(f, "fasta"):
                re_hit = _RE_RECORD_CLIP.match(record.description)
                if re_hit is None:
                    raise Exception(f"Could not parse record description: {record.description}")
                seq_id = int(re_hit.group(1))
                seq_from = int(re_hit.group(2))
                seq_to = int(re_hit.group(3))
                if (seq_from, seq_to) in out[record.id]:
                    raise Exception('Duplicate record')
                out[record.id][(seq_from, seq_to)] = str(record.seq)
        return out

    @staticmethod
    def read_file(path: Path) -> List[GraphAlignerResult]:

        out = list()
        with path.open() as f:
            for line in f.readlines():
                cols = line.strip().split('\t')

                seq_name = cols[0]
                seq_len = int(cols[1])
                seq_start = int(cols[2])
                seq_end = int(cols[3])
                seq_strand = cols[4]
                seq_path = cols[5]
                path_length = int(cols[6])
                path_start = int(cols[7])
                path_end = int(cols[8])
                n_res_matches = int(cols[9])
                aln_block_len = int(cols[10])
                map_quality = float(cols[11])

                custom_cols = dict()
                for cur_col in cols[12:]:
                    col_id, col_dtype, col_val = cur_col.split(':')
                    if col_dtype == 'i':
                        col_val = int(col_val)
                    elif col_dtype == 'f':
                        col_val = float(col_val)
                    custom_cols[col_id] = col_val

                obj = GraphAlignerResult(
                    seq_name, seq_len, seq_start, seq_end, seq_strand, seq_path,
                    path_length, path_start, path_end, n_res_matches, aln_block_len,
                    map_quality, custom_cols, None
                )
                out.append(obj)
        return out

    @staticmethod
    def version() -> str:
        """Returns the version of GraphAligner on the path."""
        proc = subprocess.Popen(
            ['GraphAligner', '--version'], encoding='utf-8',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise Exception('GraphAligner binary not found')
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
        raise Exception('GraphAligner version not found')


def _test():
    gfa = Path("/Users/aaron/git/sarand/test/spade_output/assembly_graph_with_scaffolds.gfa")
    reads = Path("/Users/aaron/git/sarand/test/aaron/amr_group_1.fasta")
    out = Path("/tmp/asd.gaf")
    clip = Path("/tmp/clip.fa")
    ga = GraphAligner.run(
        GraphAlignerParams(gfa, reads, out, preset='vg', min_align_score=1, corrected_clipped_out=clip))
    print('Done')

if __name__ == '__main__':
    _test()
