import logging
from typing import List, Tuple, Set
from collections import namedtuple

from utils.cigar import Cigar

logger = logging.getLogger()

# TODO. Do something with different letters
'''
Wildcards:
N: A, C, G or T
M: A or C
R: A or G
W: A or T
Y: C or T
S: C or G
K: G or T
H: A, C or T
B: C, G or T
V: A, C or G
D: A, G or T
'''


# TODO Only links and sequences are supported. There is nothing fancy here.
class RecordType:
    SEQUENCE = "S"
    LINK = "L"


Link: namedtuple = namedtuple("Link", ["from_name", "from_strand", "to_name", "to_strand", "cigar"])
Sequence: namedtuple = namedtuple("Sequence", ["name", "seq", "length"])


def inv_sign(sign: str) -> str:
    return ["+", "-"][sign == "+"]


def inv_link(link: Link) -> Link:
    from_name = link.to_name
    from_strand = inv_sign(link.to_strand)
    to_name = link.from_name
    to_strand = inv_sign(link.from_strand)
    cigar = str(Cigar(link.cigar).complement())
    return Link(from_name, from_strand, to_name, to_strand, cigar)


def parse_gfa(gfa_file: str) -> Tuple[List, List]:
    logger.info("Parsing file in GFA format {0}".format(gfa_file))
    sequences, links = [], []

    with open(gfa_file) as f:
        for line in f:
            record = line.strip().split()
            record_type = record[0]

            if record_type == RecordType.SEQUENCE:
                sequences.append(parse_sequence(record[1:4]))

            if record_type == RecordType.LINK:
                link = parse_link(record[1:6])
                links.append(link)

    logger.info("Done with parsing GFA file.")
    return sequences, list(links)


def parse_sequence(record: List[str]) -> Sequence:
    name: str = record[0]
    seq: str = record[1]
    length: int = int(record[2][record[2].rfind(':', 0) + 1:])
    return Sequence(name, seq, length)


def parse_link(record: List[str]) -> Link:
    from_name, from_strand, to_name, to_strand, cigar = record

    # TODO replace with proper sanity checks
    if from_name == to_name:
        print("ATTENTION, ATTENTION, ATTENTION. WE HAVE A LOOP {0} {1} {2} {3}".format(from_name, from_strand, to_name,
                                                                                       to_strand))

    return Link(from_name, from_strand, to_name, to_strand, cigar)


def write_gfa_file(gfa_fpath: str, seqs: List[Sequence], links: List[Link]) -> None:
    with open(gfa_fpath, 'w') as out:
        for seq in seqs:
            out.write("S\t{0}\t{1}\tLN:i:{2}\n".format(seq.name, seq.seq, seq.length))

        for link in links:
            out.write("L\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(link.from_name, link.from_strand,
                                                            link.to_name, link.to_strand, link.cigar))


def remove_duplicated_links(links: List[Link]) -> List[Link]:
    used = set()
    new_links: List[Link] = []

    for link in links:
        if link[0:4] not in used:
            new_links.append(link)
            used.add(link[0:4])
            used.add(inv_link(link)[0:4])

    return new_links


def save_fasta_from_graph(fasta_fpath: str, graph, vertices: Set[int]) -> None:
    logger.info("Creating fasta files " + fasta_fpath + ".")

    with open(fasta_fpath, 'w') as fasta_out:
        for v in vertices:
            fasta_out.write('>{0}\n'.format(graph.get_node_name(v)))
            length = graph.get_node_length(v)
            seq = graph.get_node_seq(v)
            for i in range(0, length, 80):
                fasta_out.write(seq[i: i + 80] + '\n')

    logger.info("Done with creating fasta file.")