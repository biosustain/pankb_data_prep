import pandas as pd
import json
import argparse
from bs4 import BeautifulSoup
from bs4.element import Tag, NavigableString, Comment
from urllib.request import urlopen
from urllib.parse import urlparse, parse_qs

IMODULONDB_BASE_URL = "https://imodulondb.org/"

def initialize_parser(parser):
    parser.description = "Parse iModulonDB for species specific info."
    parser.add_argument(
        "--im_init",
        type=str,
        required=True,
        help="iM_init generated iModulonDB json file.",
    )
    parser.add_argument(
        "--species", "-s",
        type=str,
        required=True,
        help="Species to analyze.",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output file.",
    )


def im_species(
    im_init_path,
    species,
    output_path,
):
    with open(im_init_path, "r") as f:
        org_d = json.load(f)
    
    species_name = species.replace("_", " ").replace(" E ", " ")
    if not species_name in org_d:
        with open(output_path, "w") as f:
            pass
    
    general_info = org_d[species_name]
    organism_shorthand = general_info['short']

    # gene_dl_href = f"{IMODULONDB_BASE_URL}organisms/{organism_shorthand}/annotation/gene_files.zip"
    gene_info_dl = f"{IMODULONDB_BASE_URL}organisms/{organism_shorthand}/annotation/gene_files/gene_info.csv"

    # for dataset_name, dataset_url in general_info["datasets"]:
    #     prefix = f"{IMODULONDB_BASE_URL}organisms/{general_info['short']}/{dataset_name}/"

    #     gene_info_dl = f"https://imodulondb.org/organisms/s_enterica/annotation/gene_files/gene_info.csv"
    #     # data_dl_href = f"{prefix}data_files.zip"
    #     # meta_dl_href = f"{prefix}data_files/sample_table.csv"
    #     # imtb_dl_href = f"{prefix}data_files/iM_table.csv"

    #     # print(data_dl_href)
    #     # print(gene_dl_href)
    #     # print(meta_dl_href)
    #     # print(imtb_dl_href)

    # with open(output_path, "w") as f:
    #     json.dump(org_d, f)

def run(args):
    im_species(
        args.im_init,
        args.species,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
