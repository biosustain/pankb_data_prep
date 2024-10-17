import argparse
import datetime
import logging
import os

from .nova import (
    organism,
    pangene,
    gene,
    genome,
)

logging.basicConfig(
    # filename=log_filename,
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def ask_select_mode(args):
    logging.error("Please select a mode, see --help for more info.")


def main():
    logging.info("Application started")
    parser = argparse.ArgumentParser(description=("PanKB data preparation."))
    parser.set_defaults(func=ask_select_mode)
    subparsers = parser.add_subparsers()

    modes = {
        "organism": organism,
        "pangene": pangene,
        "gene": gene,
        "genome": genome,
    }
    parsers = {}
    for x, f in modes.items():
        parsers[x] = subparsers.add_parser(x)
        f.initialize_parser(parsers[x])
        parsers[x].set_defaults(func=f.run)

    # parser.add_argument(
    #     "--version", action="version", version="%(prog)s " + __version__
    # )

    args = parser.parse_args()
    args.func(args)
