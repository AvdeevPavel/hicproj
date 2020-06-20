import argparse
import logging
import os

from __version__ import __version__
from hic_pairs.io_handler import BAMParser

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
    parser = argparse.ArgumentParser(description='Script for converting bams to .pairs',
                                     usage='pairs.py -g <graph.gfa> -b <align.bam> -o <output_dir>')

    parser.add_argument("-g", "--graph", dest="graph",
                        default=None, required=True, metavar="PATH",
                        help="Assembly graph in gfa format")

    parser.add_argument("-b", "--bam", dest="bam",
                        default=None, required=True, metavar="PATH",
                        help="Bam file with alignment of pair-end Hi-C reads")

    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="PATH", help="Output directory")

    parser.add_argument("-v", "--version", action="version", version='%(prog)s {}'.format(__version__))

    return parser.parse_args()


def main():
    # args = parse_arguments()

    # output_dir = os.path.abspath(args.out_dir)
    # if not os.path.isdir(output_dir):
    #     os.mkdir(output_dir)
    #
    # # Switching on logging
    # log_dir = os.path.join(output_dir, "logs")
    # if not os.path.isdir(log_dir):
    #     os.mkdir(log_dir)
    # log_file = os.path.join(log_dir, "pairs.log")
    # enable_logging(log_file, True, logging.DEBUG)
    enable_logging("current.log", True, logging.DEBUG)

    #
    # if not os.path.isfile(args.graph) or not os.path.isfile(args.bam):
    #     logger.error("One of the input files is unavaliable. Please, check existance of files")
    #     exit(1)
    #
    graph_file = "sb50.gfa"
    # bam_file = "two_cont.sorted.bam"
    bam_file = "algn0_mapped.bam"

    #test()
    # with open("two_cont.sam", 'r') as out:
    #
    #     # instream = iter(out)
    #     # while line is not None:
    #     #     line = next(instream, None)
    #     #     read_id = line.split('\t', 1)[0] if line else None
    #     #     if read_id != '@HD':
    #     #         break
    #
    #     streaming_classify(out)

    parser = BAMParser(bam_file, graph_file)
    parser.parse()
    parser.write_statistics_info_into_file("stat1.txt")


if __name__ == "__main__":
    main()
