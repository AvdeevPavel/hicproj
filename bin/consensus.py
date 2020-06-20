from typing import List, Dict

from gfa_graph.asm_graph import AsmGraph, build_graph
from gfa_graph.io_handler import parse_gfa, remove_duplicated_links
from gfa_graph.superbubbles import detect_all_superbubbles, hide_complex_bubbles, get_complex_bubbles, \
    get_simple_bubbles

interesting_tags = {'tp',  # A	Type of aln: P/primary, S/secondary and I,i/inversion
                    'NM',  # i	Total number of mismatches and gaps in the alignment
                    'cs',  # Z	Difference string
                    'dv'
                    }


class PAFRecord(object):
    def __init__(self, inpt_str: str):
        elems = inpt_str.strip().split('\t')
        self.info = dict()
        self.info["qname"] = elems[0]  # 1 str query sequence name
        self.info["qlength"] = elems[1]  # 2	int	Query sequence length
        self.info["qstart"] = elems[2]  # 3	int	Query start coordinate (0-based)
        self.info["qend"] = elems[3]  # 4  int	Query end coordinate (0-based)
        self.info["same_strand"] = elems[4]  # 5	char	‘+’ if query/target on the same strand; ‘-’ if opposite
        self.info['tname'] = elems[5]  # 6	string	Target sequence name
        self.info['tlength'] = elems[6]  # 7	int	Target sequence length
        self.info["tstart"] = elems[7]  # 8	int	Target start coordinate on the original strand
        self.info["tend"] = elems[8]  # 9	int	Target end coordinate on the original strand
        self.info["nmatch"] = elems[9]  # 10	int	Number of matching bases in the mapping
        self.info["nbases"] = elems[10]  # 11	int	Number bases, including gaps, in the mapping
        self.info["mapq"] = elems[11]  # 12	int	Mapping quality (0-255 with 255 for missing)
        self._parse_tags(elems[12:])

    def _parse_tags(self, tags: List):
        self.tags = dict()
        for tag in tags:
            tag_name = tag[0:tag.find(':')]
            tag_type = tag[tag.find(':') + 1:tag.rfind(':')]
            tag_val = tag[tag.rfind(':') + 1:]
            # TODO add converter depending on tag type
            self.tags[tag_name] = (tag_type, tag_val)

    def has_tag(self, name):
        return name in self.tags

    def get_tag_value(self, name):
        return self.tags[name][1]

    def _tags_to_string(self) -> str:
        result: str = ''
        for tag_name, tag_prop in self.tags.items():
            tag_type, tag_val = tag_prop
            result += "{0}:{1}:{2}\t".format(tag_name, tag_type, tag_val)
        return result

    def to_string(self) -> str:
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}".format(
            self.info["qname"],
            self.info["qlength"],
            self.info["qstart"],
            self.info["qend"],
            self.info["same_strand"],
            self.info['tname'],
            self.info['tlength'],
            self.info["tstart"],
            self.info["tend"],
            self.info["nmatch"],
            self.info["nbases"],
            self.info["mapq"],
            self._tags_to_string()
        )


class PAFParser(object):
    def __init__(self, paf_file: str) -> None:
        self.paf_file = paf_file
        self.total_counter = 0
        self.primary_counter = 0
        self.secondary_counter = 0
        self.inversion_counter = 0

        self.expected_counter = 0
        self.unexpected_counter = 0

    def parse(self, bubbly_info: Dict) -> None:
        with open(self.paf_file, 'r') as inp:
            for line in inp:
                self.total_counter += 1
                record = PAFRecord(line)

                if record.has_tag('tp'):
                    if record.get_tag_value('tp') == 'P':
                        self.primary_counter += 1
                    elif record.get_tag_value('tp') == 'S':
                        self.secondary_counter += 1
                    else:
                        self.inversion_counter += 1

                if bubbly_info[record.info["qname"]] == bubbly_info[record.info["tname"]]:
                    if record.get_tag_value('tp') == 'P':
                        self.expected_counter += 1
                    else:
                        self.unexpected_counter += 1

    def write_stat(self, out_file: str) -> None:
        with open(out_file, 'w') as out:
            out.write("Total number of records = {0}\n".format(self.total_counter))
            out.write("Total number of primary alignments = {0}\n".format(self.primary_counter))
            out.write("Total number of secondary alignments = {0}\n".format(self.secondary_counter))
            out.write("Total number of inversion-like alignments = {0}\n".format(self.inversion_counter))
            out.write("\n")
            out.write("Total number of primary and in bubble = {0}\n".format(self.expected_counter))
            out.write("Total number of secondary and in bubble = {0}\n".format(self.unexpected_counter))


def get_bubbly_map(graph_path: str) -> Dict:
    seqs, links = parse_gfa(graph_path)
    links = remove_duplicated_links(links)

    graph: AsmGraph = build_graph(seqs, links)
    superbubbles = detect_all_superbubbles(graph)
    hide_complex_bubbles(graph, get_complex_bubbles(superbubbles))

    result = dict()
    for bubble, ind in zip(get_simple_bubbles(superbubbles), range(len(get_simple_bubbles(superbubbles)))):
        b_nodes = list(bubble.internal_nodes())
        assert len(b_nodes) == 2
        result[graph.get_node_name(b_nodes[0].id)] = ind
        result[graph.get_node_name(b_nodes[1].id)] = ind
    return result


def main():
    paf_parser = PAFParser("../algn1.paf")
    paf_parser.parse(get_bubbly_map("sb50.gfa"))
    paf_parser.write_stat("paf_stat.txt")


if __name__ == '__main__':
    main()
