import argparse
import datetime
import logging
import os

from . import (
    isolation_source,
    eggnog_summary,
    all_summary,
    species_summary,
    family_summary,
    full_summary,
    alleleome_dimensions,
    mash_list,
    heatmap,
    cog_data,
    cog_distribution,
    heaps_law,
    gene_locustag,
    genome_page,
    landing_page,
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
        "isosource": isolation_source,
        "eggnog": eggnog_summary,
        "full": full_summary,
        "family": family_summary,
        "species": species_summary,
        "all": all_summary,
        "alleleome_dimensions": alleleome_dimensions,
        "mash": mash_list,
        "heatmap": heatmap,
        "cog": cog_data,
        "distribution": cog_distribution,
        "heaps": heaps_law,
        "locustag": gene_locustag,
        "genome": genome_page,
        "landing": landing_page,
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
