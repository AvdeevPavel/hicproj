import logging
import sys
from collections import namedtuple
from typing import Dict, List, Set

from gfa_graph.asm_graph import AsmGraph
from gfa_graph.graph_structs import DirectedNode, Arc, possible_orientations

logger = logging.getLogger()

Range: namedtuple = namedtuple("Range", ["start", "end"])


# TODO add a lot of tests A LOT
class Superbubble(object):
    """
        Miniasm description -- https://academic.oup.com/bioinformatics/article/32/14/2103/1742895
    """

    def __init__(self, sv: DirectedNode, max_length: int = -1, max_cnt: int = -1):
        self._start_vertex: DirectedNode = sv
        self._end_vertex: DirectedNode = sv
        self._n_vertices: int = 0

        # stores longest path
        self._superbubble_vertices: Dict[DirectedNode, int] = dict()
        # stores previous optimal path
        self._heaviest_backtrace: Dict[DirectedNode, DirectedNode] = dict()

        self._longest_length: int = sys.maxsize
        self._max_cnt = sys.maxsize

    def start_vertex(self) -> DirectedNode:
        return self._start_vertex

    def end_vertex(self) -> DirectedNode:
        return self._end_vertex

    def update_length(self, graph: AsmGraph, v: DirectedNode, arc: Arc) -> None:
        node_length = graph.get_node_length(arc.dn.id)
        overlap_len = arc._cigar.size_on_query()

        assert overlap_len < node_length

        if arc.dn not in self._superbubble_vertices or \
                self._superbubble_vertices[arc.dn] < self._superbubble_vertices[v] + node_length - overlap_len:
            self._superbubble_vertices[arc.dn] = self._superbubble_vertices[v] + node_length - overlap_len
            self._heaviest_backtrace[arc.dn] = v

    def nodes(self):
        return self._superbubble_vertices.keys()

    # TODO Problem with tips in bubbles and if end_vertex in itself is tip. (I patched but need to be tested)
    def detect(self, graph) -> bool:
        if graph.outgoing_degree(self._start_vertex) < 2:
            logger.debug("Vertex {0} has outgoing degree < 2. No bubbles starting from here".format(
                graph.get_dir_node_name(self._start_vertex)))
            return False

        self._superbubble_vertices[self._start_vertex] = graph.get_node_length(self._start_vertex.id)
        self._heaviest_backtrace[self._start_vertex] = DirectedNode()
        self._n_vertices += 1

        unvis_incoming_edges: Dict[DirectedNode, int] = dict()
        stack: List = [self._start_vertex]
        p: int = 0

        while len(stack) != 0:
            v = stack.pop()
            logger.debug("Working with {0} vertex in bubble".format(graph.get_dir_node_name(v)))

            if self._n_vertices + 1 > self._max_cnt:
                return False
            self._n_vertices += 1

            # logger.debug("Outgoing vertices are {0}".format([graph.get_node_name(l._dn.id)
            #                                                  for l in graph.outgoing(v)]))
            for l in graph.outgoing(v):
                neighbour_v: DirectedNode = l.dn

                if neighbour_v == self._start_vertex:
                    logger.debug("We have a cycle that includes start vertex")
                    return False

                if neighbour_v not in self._superbubble_vertices:
                    logger.debug("Registering new vertex for consideration {0}".format(
                        graph.get_dir_node_name(neighbour_v)))
                    unvis_incoming_edges[neighbour_v] = graph.incoming_degree(neighbour_v)
                    p += 1

                self.update_length(graph, v, l)

                unvis_incoming_edges[neighbour_v] -= 1

                if unvis_incoming_edges[neighbour_v] == 0:
                    p -= 1
                    if (p == 0 and graph.outgoing_degree(neighbour_v) == 0) or \
                            graph.outgoing_degree(neighbour_v) != 0:
                        stack.append(neighbour_v)

            # logger.debug("Stack length {0} and #unvisited vertices {1}".format(len(stack), p))
            # logger.debug("Registering new vertex for consideration {0}".format(

            if len(stack) == 1 and p == 0:
                self._end_vertex = stack.pop()
                logger.debug("Sink was found {0}".format(graph.get_dir_node_name(self._end_vertex)))

        if p == 0:
            # logger.debug("Sink was found1 {0}".format(self._graph.get_dir_node_name(self._end_vertex)))
            # logger.debug("Outgoing vertices are {0}".format([self._graph.get_dir_node_name(l._dn)
            #                                                  for l in self._graph.incoming(self._end_vertex)]))

            max_l: int = 0
            parent: DirectedNode = DirectedNode()
            for l in graph.incoming(self._end_vertex):
                assert l.dn in self._superbubble_vertices
                if max_l < self._superbubble_vertices[l.dn]:
                    max_l = self._superbubble_vertices[l.dn]
                    parent: DirectedNode = l.dn

            assert parent != DirectedNode()
            self._superbubble_vertices[self._end_vertex] = max_l
            self._heaviest_backtrace[self._end_vertex] = parent
            return True
        else:
            return False

    def is_simple_bubble(self) -> int:
        return len(self._superbubble_vertices) == 4

    def get_nodes_in_longest_path(self) -> Set[DirectedNode]:
        assert self._end_vertex != self._start_vertex
        result: Set[DirectedNode] = set()
        v: DirectedNode = self._end_vertex
        while v != self._start_vertex:
            logger.debug("Added to longest path a node {0}".format(v.id))
            result.add(v)
            v = self._heaviest_backtrace[v]
        result.add(self._start_vertex)
        return result


def detect_all_superbubbles(graph: AsmGraph) -> List[Superbubble]:
    logger.info("Finding bubbles in graph")
    in_bubble: Set[DirectedNode] = set()
    detected_bubbles: Dict[DirectedNode, Superbubble] = dict()

    for id in graph.node_ids():
        for o in possible_orientations():
            v: DirectedNode = DirectedNode(id, o)
            logger.debug("Looking at directed node {0}".format(graph.get_dir_node_name(v)))

            if v in in_bubble:
                logger.debug("Not considering. Was part of bubble.".format(id))
                continue

            finder: Superbubble = Superbubble(v)
            if finder.detect(graph):
                logger.info(
                    "Found superbubble between {0} and {1}".format(graph.get_dir_node_name(finder.start_vertex()),
                                                                   graph.get_dir_node_name(finder.end_vertex())))

                detected_bubbles[finder.start_vertex()] = finder

                for dn in finder.nodes():
                    if dn != finder.end_vertex() and dn != finder.start_vertex():
                        if dn in detected_bubbles:
                            detected_bubbles.pop(dn)
                        elif dn.complement() in detected_bubbles:
                            detected_bubbles.pop(dn.complement())

                    if dn != finder.end_vertex():
                        in_bubble.add(dn)

                    if dn != finder.start_vertex():
                        in_bubble.add(dn.complement())

    logger.info("Done with finding bubbles in graph")
    return list(detected_bubbles.values())


def get_complex_bubbles(bubbles: List[Superbubble]) -> List[Superbubble]:
    return [x for x in bubbles if not x.is_simple_bubble()]


def get_simple_bubbles(bubbles: List[Superbubble]) -> List[Superbubble]:
    return [x for x in bubbles if x.is_simple_bubble()]


def hide_complex_bubbles(graph: AsmGraph, complex_bubbles: List[Superbubble]):
    logger.info("Deleting nodes in superbubbles that are not in the longest path")
    for sb in complex_bubbles:
        path_nodes = set(sb.get_nodes_in_longest_path())
        for dn in sb.nodes():
            if dn not in path_nodes:
                graph.delete_node(dn)
    logger.info("Done with deleting nodes in superbubbles that are not in the longest path")
