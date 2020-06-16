import logging
from collections import defaultdict
from typing import List, Dict, Set
from gfa_graph.graph_structs import NodeInfo, Arc, DirectedNode, NodeOrientation, string2orientation, \
    possible_orientations, orientation2string
from gfa_graph.io_handler import Sequence, Link, inv_link
from utils.cigar import Cigar

logger = logging.getLogger()


class AsmGraph(object):
    def __init__(self):
        self._nodes: Dict[int, NodeInfo] = dict()
        self._node2id: Dict[str, int] = dict()
        self._max_node_id: int = 0

        self._outgoing_arcs: Dict[int, List[Arc]] = defaultdict(list)
        self._incoming_arcs: Dict[int, List[Arc]] = defaultdict(list)

    def add_node(self, name, seq, length) -> None:
        self._nodes[self._max_node_id] = NodeInfo(name=name, seq=seq, length=length)
        self._node2id[name] = self._max_node_id
        self._max_node_id += 1

    def get_node_name(self, id: int) -> str:
        assert id in self._nodes
        return self._nodes[id].name

    def get_node_seq(self, id: int) -> str:
        assert id in self._nodes
        return self._nodes[id].seq

    def get_node_length(self, id: int) -> int:
        assert id in self._nodes
        return self._nodes[id].length

    def get_dir_node_name(self, v: DirectedNode) -> str:
        assert v.id in self._nodes
        return orientation2string(v.o) + self._nodes[v.id].name

    def node_ids(self):
        return self._nodes.keys()

    def delete_node(self, dv: DirectedNode):
        name = self.get_node_name(dv.id)
        logger.debug("Deleting node {0} with id {1}".format(name, dv.id))

        self._node2id.pop(name)
        self._nodes.pop(dv.id)

        for id in self.adjacent_nodes(dv.id):
            if id in self._incoming_arcs:
                self._incoming_arcs[id] = list(filter(lambda x: x.dn.id != dv.id, self._incoming_arcs[id]))

            if id in self._outgoing_arcs:
                self._outgoing_arcs[id] = list(filter(lambda x: x.dn.id != dv.id, self._outgoing_arcs[id]))

        if dv.id in self._incoming_arcs:
            self._incoming_arcs.pop(dv.id)

        if dv.id in self._outgoing_arcs:
            self._outgoing_arcs.pop(dv.id)

        logger.debug("Done with deleting node {0} with id {1}".format(name, dv.id))

    def get_node_id(self, name: str) -> int:
        assert name in self._node2id
        return self._node2id[name]

    def add_link(self, source_id: int, source_orient: NodeOrientation,
                 sink_id: int, sink_orient: NodeOrientation, cigar: Cigar) -> None:
        if source_orient == NodeOrientation.FORWARD:
            if sink_orient == NodeOrientation.FORWARD:
                self._outgoing_arcs[source_id].append(Arc(DirectedNode(sink_id, NodeOrientation.FORWARD), cigar))
                self._incoming_arcs[sink_id].append(Arc(DirectedNode(source_id, NodeOrientation.FORWARD), cigar))
            else:
                self._outgoing_arcs[source_id].append(Arc(DirectedNode(sink_id, NodeOrientation.REVERSE), cigar))
                self._outgoing_arcs[sink_id].append(Arc(DirectedNode(source_id, NodeOrientation.REVERSE),
                                                        cigar.complement()))
        else:
            if sink_orient == NodeOrientation.FORWARD:
                self._incoming_arcs[source_id].append(Arc(DirectedNode(sink_id, NodeOrientation.REVERSE),
                                                          cigar.complement()))
                self._incoming_arcs[sink_id].append(Arc(DirectedNode(source_id, NodeOrientation.REVERSE), cigar))
            else:
                self._incoming_arcs[source_id].append(Arc(DirectedNode(sink_id, NodeOrientation.FORWARD),
                                                          cigar.complement()))
                self._outgoing_arcs[sink_id].append(Arc(DirectedNode(source_id, NodeOrientation.FORWARD),
                                                        cigar.complement()))

    def number_of_nodes(self):
        return len(self._nodes)

    def incoming(self, dn: DirectedNode) -> List[Arc]:
        return self._incoming_arcs[dn.id] if dn.o == NodeOrientation.FORWARD \
            else complement_links(self._outgoing_arcs[dn.id])

    def incoming_degree(self, dn: DirectedNode) -> int:
        return len(self._incoming_arcs[dn.id]) if dn.o == NodeOrientation.FORWARD \
            else len(self._outgoing_arcs[dn.id])

    def outgoing(self, dn: DirectedNode) -> List[Arc]:
        return self._outgoing_arcs[dn.id] if dn.o == NodeOrientation.FORWARD \
            else complement_links(self._incoming_arcs[dn.id])

    def outgoing_degree(self, dn: DirectedNode) -> int:
        return len(self._outgoing_arcs[dn.id]) if dn.o == NodeOrientation.FORWARD \
            else len(complement_links(self._incoming_arcs[dn.id]))

    def adjacent_nodes(self, id: int):
        assert id < self._max_node_id
        return [l.dn.id for l in self._incoming_arcs[id]] + [l.dn.id for l in self._outgoing_arcs[id]]

    def create_gfa_info(self):
        logger.info("Converting assembly graph to GFA records")
        sequences: List[Sequence] = [node.to_sequence_record() for node in self._nodes.values()]

        def arc_to_link(sdn: DirectedNode, arc: Arc) -> Link:
            return Link(self._nodes[sdn.id].name, orientation2string(sdn.o),
                        self._nodes[arc.dn.id].name, orientation2string(arc.dn.o), str(arc._cigar))

        # TODO need to rethink it and get rid of considered links
        links: List[Link] = []
        considered_links: Set[Link] = set()
        for id, n in self._nodes.items():
            logger.debug("Outputting links for node {0} with {1}".format(n.name, id))
            for o in possible_orientations():
                v: DirectedNode = DirectedNode(id, o)

                for arc in self.outgoing(v):

                    slink = arc_to_link(v, arc)
                    if slink[0:4] not in considered_links:
                        links.append(slink)
                        considered_links.add(slink[0:4])
                        considered_links.add(inv_link(slink)[0:4])

        logger.info("Done with converting assembly graph to GFA records")
        return sequences, links


def complement_links(links: List[Arc]) -> List[Arc]:
    return [l.complement() for l in links]


def build_graph(sequences: List[Sequence], links: List[Link]) -> AsmGraph:
    logger.info("Building assembly graph")
    graph = AsmGraph()

    logger.info("Adding nodes to graph")
    for seq_r in sequences:
        graph.add_node(seq_r.name, seq_r.seq, seq_r.length)

    logger.info("Adding links to graph")
    for link in links:
        graph.add_link(graph.get_node_id(link.from_name), string2orientation(link.from_strand),
                       graph.get_node_id(link.to_name), string2orientation(link.to_strand), Cigar(link.cigar))

    logger.info("Done with building assembly graph")
    return graph
