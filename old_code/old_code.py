import os
import re
import sys
from typing import List, Dict

import networkx as nx
import yaml


def prepare_additional_files(gfa_fpath: str, fasta_fpath: str, contsz_fpath: str, coord_fpath: str,
                             res_enzymes: List[str]) -> None:
    """
    TODO
    :param gfa_fpath:
    :param fasta_fpath:
    :param res_enzymes:
    :param coord_fpath:
    :param contsz_fpath:
    :return:
    """
    print("Creating fasta file from gfa file " + gfa_fpath + " ...")

    # TODO check existing of the file (catch exception)
    with open(gfa_fpath, 'r') as gfa_in, \
            open(fasta_fpath, 'w') as fasta_out, \
            open(contsz_fpath, 'w') as contsz_out, \
            open(coord_fpath, 'w') as coord_out:
        for line in gfa_in:
            record_type = line[0]
            if record_type == 'S':
                # parse vertex of gfa (segment)
                fs = line.split()
                name, seq, slen = fs[1], fs[2], fs[3]
                length = int(slen[slen.rfind(':', 0) + 1:])
                assert length == len(seq)

                # write unitig to fasta
                fasta_out.write('> {0}\n'.format(name))
                for i in range(0, length, 80):
                    fasta_out.write(seq[i: i + 80] + '\n')

                # write contig sizes
                contsz_out.write('{0}\t{1}\n'.format(name, length))

                # TODO what if seq contains N's or other wildcards in Seq (in enzymes we can deal before)?
                # TODO What is about reverse complement?
                # TODO It can be optimized with proper implementation of aho-corasick algo
                if res_enzymes[0] != 'DNASE':

                    if len(res_enzymes) == 1:
                        regex = re.compile(res_enzymes[0], re.IGNORECASE)
                    else:
                        regex = re.compile('|'.join(res_enzymes), re.IGNORECASE)

                    poss = [str(m.start()) for m in regex.finditer(seq)]
                    coord_out.write('{0}\t{1}\n'.format(name, '\t'.join(poss)))

def parse_file_with_hic_libs(yaml_file: str) -> Dict:
    """
    TODO
    :param yaml_file:
    :return:
    """

    # TODO check existing of the file
    with open(yaml_file, 'r') as inp:
        libs = yaml.full_load(inp)

        for lib in libs:
            final_enzymes: List = []

            for x in lib['res_sites']:
                if x == "DNASE":
                    final_enzymes: List = ["DNASE"]
                    break
                elif not re.match("^[ACGTN]+$", x):
                    # TODO redo on exceptions
                    print("Error, enzyme should be restriction site sequence (e.g. AACTT) "
                          "not enzyme name or DNASE, you input %s".format(x), sys.stderr)
                    sys.exit(1)
                else:
                    if 'N' in x:
                        final_enzymes.append(x.replace('N', 'G'))
                        final_enzymes.append(x.replace('N', 'A'))
                        final_enzymes.append(x.replace('N', 'T'))
                        final_enzymes.append(x.replace('N', 'C'))
                    else:
                        final_enzymes.append(x)

            lib['res_sites'] = final_enzymes

        return libs


def generate_files_for_snakemake(libs: Dict, output_dir: str) -> None:
    """
    TODO
    :param libs:
    :return:
    """
    for lib in libs:
        fdir = os.path.join(output_dir, lib["name"])

        if not os.path.isdir(fdir):
            os.mkdir(fdir)

        with open(os.path.join(fdir, "snake_desc.bed"), 'w') as out_snake:
            for read_pairs, ind in zip(lib["fastqs"], range(0, len(lib["fastqs"]))):
                out_snake.write("{0} {1} {2}")


def form_fasta_file(gfa_fpath: str, fasta_fpath: str, contsz_fpath: str) -> None:
    """
    TODO
    :param gfa_fpath:
    :param fasta_fpath:
    :return:
    """
    print("Creating fasta file from gfa file " + gfa_fpath + " ...")

    # TODO check existing of the file (catch exception)
    with open(gfa_fpath, 'r') as gfa_in, \
            open(fasta_fpath, 'w') as fasta_out, \
            open(contsz_fpath, 'w') as contsz_out:
        for line in gfa_in:
            record_type = line[0]
            if record_type == 'S':
                # parse vertex of gfa (segment)
                fs = line.split()
                name, seq, slen = fs[1], fs[2], fs[3]
                length = int(slen[slen.rfind(':', 0) + 1:])

                assert length == len(seq)

                # write unitig to fasta
                fasta_out.write('> {0}\n'.format(name))
                for i in range(0, length, 80):
                    fasta_out.write(seq[i: i + 80] + '\n')

                # write contig sizes
                contsz_out.write('{0}\t{1}\n'.format(name, length))


def statistics_about_gfa(gfa_fpath: str, stat_fpath: str) -> None:
    print("Parsing " + gfa_fpath + "...")

    tot_length = 0
    # graph = nx.Graph()
    c_s, c_l, c_a = 0, 0, 0
    with open(gfa_fpath) as f:
        for line in f:
            record_type = line[0]
            if record_type == 'S':
                fs = line.split()
                name, seq, len = fs[1], fs[2], int(fs[3].split(":")[2])
                # graph.add_node(name, name=name, length=len)
                # graph.add_node(name + "$REV", name=name, length=len)
                # graph.add_node(name + "$FOW", name=name, length=len)
                tot_length += len
                c_s += 1
            elif record_type == 'L':
                c_l += 1
            elif record_type == 'A':
                # IGNORE FOR NOW
                c_a += 1
            elif record_type == 'E':
                # IGNORE FOR NOW
                pass

    with open(stat_fpath, 'w') as sout:
        sout.write("Total segment length {0}\n".format(tot_length))
        sout.write("Average segment length {0}\n".format(float(tot_length) / c_s))
        sout.write("Number of segments {0}\n".format(c_s))
        sout.write("Number of edges {0}\n".format(c_l))
        sout.write("Number of arcs(?) {0}\n".format(c_a))


# def parse_gfa(gfa_fpath: str): # , out_path: str):
#     print("Parsing " + gfa_fpath + "...")
#     with open(gfa_fpath, 'r') as f: #, open(out_path, 'w') as out:
#         for line in f:
#             record_type = line[0]
#             if record_type == 'S':
#                 pass
#                 # out.write(line)
#                 # fs = line.split()
#                 # name, seq, len = fs[1], fs[2], int(fs[3].split(":")[2])
#             elif record_type == 'L':
#                 pass
#                 # out.write(line)
#                 # pass
#             elif record_type == 'A':
#                 fs = line.split()
#                 name, l, strand, rname, x, y, n, tag = fs[1], fs[2], fs[3], fs[4], fs[5], fs[6], fs[7], fs[8].split(":")[2]
#                 if name == "utg027119l": # "utg037337l":
#                     print(line)


def get_induced_subgraph(gfa_fpath: str):
    print("Parsing " + gfa_fpath + "...")

    graph = nx.Graph()

    with open(gfa_fpath) as f:
        for line in f:
            record_type = line[0]
            if record_type == 'S':
                fs = line.split()
                name, seq, length = fs[1], fs[2], int(fs[3].split(":")[2])
                # graph.add_node(name, name=name, length=len)
                graph.add_node(name + "$REV", name=name, length=length)
                graph.add_node(name + "$FOW", name=name, length=length)
                # tot_length += int(slen[slen.rfind(':', 0) + 1:])
                # c_s += 1
            elif record_type == 'L':
                # TODO cigar only with M right now. Support more letters
                _, from_name, from_orient, to_name, to_orient, cigar = line.split()[:6]

                if from_orient == "+" and to_orient == "+":
                    graph.add_edge(from_name + "$FOW", to_name + "$FOW", weight=1)
                    graph.add_edge(to_name + "$REV", from_name + "$REV", weight=1)
                if from_orient == "+" and to_orient == "-":
                    graph.add_edge(from_name + "$FOW", to_name + "$REV", weight=1)
                    graph.add_edge(to_name + "$FOW", from_name + "$REV", weight=1)
                if from_orient == "-" and to_orient == "+":
                    graph.add_edge(from_name + "$REV", to_name + "$FOW", weight=1)
                    graph.add_edge(to_name + "$REV", from_name + "$FOW", weight=1)
                if from_orient == "-" and to_orient == "-":
                    graph.add_edge(from_name + "$REV", to_name + "$REV", weight=1)
                    graph.add_edge(to_name + "$FOW", from_name + "$FOW", weight=1)

    print("The number of connected componets {0}".format(nx.number_connected_components(graph)))
    min_comp = list(filter(lambda x: len(x) > 200, [c for c in sorted(nx.connected_components(graph), key=len)]))
    # min_comp = min(nx.connected_components(graph), key=len)
    print("The number of connected componets {0}".format(len(min_comp)))
    print("Minimum size component is {0}".format([len(x) for x in min_comp[:10]]))

    print("Save subgraphs")
    vertex_set = set()
    for n1 in min_comp[0]:
        vertex_set.add(n1.split('$')[0])

        # vertex_set.add(e2.split('$')[0])
    # for e1, e2 in nx.dfs_edges(graph, source=list(graph.nodes)[0], depth_limit=min_depth):
    #     vertex_set.add(e1.split('$')[0])
    #     vertex_set.add(e2.split('$')[0])

    with open(gfa_fpath, 'r') as f, open("subgraph.gfa", 'w') as out:
        for line in f:
            record_type = line[0]
            if record_type == 'S':
                name = line.split()[1]
                if name in vertex_set:
                    out.write(line)
            elif record_type == 'L':
                _, from_name, from_orient, to_name, to_orient, cigar = line.split()[:6]
                if from_name in vertex_set and to_name in vertex_set:
                    out.write(line)

    # if cigar.find('I') != -1 or cigar.find('')
    # overlap_operations = re.split('(\d+)', cigar)
    # print(overlap_operations)
    # elif record_type == 'A':


