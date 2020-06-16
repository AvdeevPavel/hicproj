from enum import Enum
from typing import List
from utils.cigar import Cigar
from gfa_graph.io_handler import Sequence

vertex_id_t = int


class NodeOrientation(Enum):
    FORWARD = 0
    REVERSE = 1


def possible_orientations() -> List:
    return [NodeOrientation.FORWARD, NodeOrientation.REVERSE]


def get_opposite_orientation(o: NodeOrientation) -> NodeOrientation:
    return NodeOrientation.REVERSE if o == NodeOrientation.FORWARD else NodeOrientation.FORWARD


def orientation2string(o: NodeOrientation) -> str:
    return '+' if o == NodeOrientation.FORWARD else '-'


def string2orientation(c: str) -> NodeOrientation:
    return NodeOrientation.FORWARD if c == '+' else NodeOrientation.REVERSE


class NodeInfo(object):
    def __init__(self, name: str = '', length: int = 0, seq: str = '') -> None:
        self.name: str = name
        self.length: int = length
        self.seq: str = seq

    def to_sequence_record(self) -> Sequence:
        return Sequence(self.name, self.seq, self.length)


class DirectedNode(object):
    def __init__(self, id: vertex_id_t = '', o: NodeOrientation = NodeOrientation.FORWARD) -> None:
        self.id: vertex_id_t = id
        self.o: NodeOrientation = o

    def complement(self):
        return DirectedNode(self.id, get_opposite_orientation(self.o))

    def __hash__(self):
        return 2 * self.id + (0 if self.o == NodeOrientation.FORWARD else 1)

    def __eq__(self, other: "DirectedNode"):
        return self.id == other.id and self.o == other.o

    def __ne__(self, other):
        return not self.__eq__(other)


class Arc(object):
    def __init__(self, dn: DirectedNode, cigar: Cigar):
        self.dn: DirectedNode = dn
        self._cigar: Cigar = cigar

    def complement(self):
        return Arc(self.dn.complement(), self._cigar.complement())
