from typing import List

from sarand.external.graph_aligner import GraphAlignerResult


class GraphAlignmentOutput:
    """Interface to provide a modular way to retrieve the core dataset required
    for the graph alignment step.
    """

    def __init__(
            self,
            nodes: List[int],
            orientations: List[str],
            start_pos: int,
            end_pos: int
    ):
        self.nodes: List[int] = nodes
        self.orientation: List[str] = orientations
        self.start_pos: int = start_pos
        self.end_pos: int = end_pos

    @classmethod
    def from_graph_aligner(cls, ga: List[GraphAlignerResult]):

        return
