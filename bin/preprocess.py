#!/usr/bin/env python3

import argparse
import logging
import os

from __version__ import __version__
from gfa_graph.asm_graph import build_graph, AsmGraph

from gfa_graph.io_handler import parse_gfa, write_gfa_file, save_fasta_from_graph, remove_duplicated_links
from gfa_graph.superbubbles import detect_all_superbubbles, get_complex_bubbles, hide_complex_bubbles, \
    get_simple_bubbles


logger = logging.getLogger()


def enable_logging(log_file: str, overwrite: bool, log_level=logging.INFO) -> None:
    """
    Turns on logging, sets debug levels and assigns a log file.
    Args:
            log_file (str): The path to log file
            overwrite (bool): The indicator for overwriting existing log file
            log_level:
    Returns:
         no value
    """
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(log_level)
    logger.addHandler(console_log)

    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    if overwrite:
        open(log_file, "w").close()

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(log_level)
    logger.addHandler(file_handler)


def parse_arguments():
    """
    TODO
    Parse input arguments
    :return:
    """
    parser = argparse.ArgumentParser(description='Script for detecting bubbles and forming fasta files for alignments',
                                     usage='graph.py -g <graph.gfa> -o <output_dir>')

    parser.add_argument("-g", "--graph", dest="graph",
                        default=None, required=True, metavar="PATH",
                        help="Assembly graph in gfa format")

    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="PATH", help="Output directory")

    parser.add_argument("-v", "--version", action="version", version='%(prog)s {}'.format(__version__))

    return parser.parse_args()


def main():
    args = parse_arguments()

    # Checking output dir
    output_dir = os.path.abspath(args.out_dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Switching on logging
    log_dir = os.path.join(output_dir, "logs")
    if not os.path.isdir(log_dir):
        os.mkdir(log_dir)
    log_file = os.path.join(log_dir, "preprocess.log")
    enable_logging(log_file, True, logging.DEBUG)

    # Checking existence of the graph
    graph_path = os.path.abspath(args.graph)
    if not os.path.isfile(graph_path):
        logger.error("File {0} with graph is unavaliable.".format(graph_path))
        exit(1)

    graph_out_dir = os.path.join(output_dir, "inp_graph")
    if not os.path.isdir(graph_out_dir):
        os.mkdir(graph_out_dir)

    seqs, links = parse_gfa(graph_path)
    links = remove_duplicated_links(links)
    write_gfa_file(os.path.join(graph_out_dir, "graph.gfa"), seqs, links)

    # Finding bubbles
    graph: AsmGraph = build_graph(seqs, links)
    superbubbles = detect_all_superbubbles(graph)
    hide_complex_bubbles(graph, get_complex_bubbles(superbubbles))

    # Saving two fastas (bubble and non-bubble contigs)
    # TODO deal with overlaps
    fasta_out_dir = os.path.join(output_dir, "fastas")
    if not os.path.isdir(fasta_out_dir):
        os.mkdir(fasta_out_dir)

    bubbly_out_dir = os.path.join(fasta_out_dir, "bubbly_utg")
    if not os.path.isdir(bubbly_out_dir):
        os.mkdir(bubbly_out_dir)

    non_bubbly_out_dir = os.path.join(fasta_out_dir, "non_bubbly_utg")
    if not os.path.isdir(non_bubbly_out_dir):
        os.mkdir(non_bubbly_out_dir)

    utgs = {n.id for bubble in get_simple_bubbles(superbubbles) for n in bubble.nodes()}
    save_fasta_from_graph(os.path.join(bubbly_out_dir, "seqs.fa"), graph, utgs)
    utgs = set(graph.node_ids()).difference(utgs)
    save_fasta_from_graph(os.path.join(non_bubbly_out_dir, "seqs.fa"), graph, utgs)

    seqs, links = graph.create_gfa_info()
    write_gfa_file(os.path.join(output_dir, "graph_wo_csb.gfa"), seqs, links)


if __name__ == "__main__":
    main()
