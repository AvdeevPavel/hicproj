import argparse
import logging
import pysam

from __version__ import __version__

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
    args = parse_arguments()



if __name__ == "__main__":
    main()
