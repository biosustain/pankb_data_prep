import argparse
import pandas as pd
import requests
from pathlib import Path
import json
import os
from .imodulon_data import IMODULON_DATA

def touch(path):
    with open(path, 'a'):
        os.utime(path, None)

def download_file(url, path):
    r = requests.get(url)
    with open(path, "wb") as f:
        f.write(r.content)

def initialize_parser(parser):
    parser.description = "Download iModulon data for species."
    parser.add_argument(
        "name",
        type=str,
        help="Family or analysis name to output data for.",
    )
    parser.add_argument(
        "--samples",
        type=str,
        required=True,
        help="Samples input csv.",
    )
    parser.add_argument(
        "--genome_list",
        type=str,
        required=True,
        help="Genome list output.",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output directory.",
    )

def organism_data(
    analysis_name,
    samples_path,
    genome_list_path,
    output_path,
):
    if not analysis_name in IMODULON_DATA:
        touch(genome_list_path)
        Path(output_path).mkdir(parents=True, exist_ok=True)
    
    df_samples = pd.read_csv(samples_path, index_col="genome_id")
    
    cur_organism = IMODULON_DATA[analysis_name]
    genomes = []
    for genome in cur_organism.keys():
        if genome in df_samples.index:
            genomes.append(genome)
        else:
            print(f"Genome {genome} present in iModulonDB, but not PanKB")
    df = pd.DataFrame({"genome_id": genomes})
    df.to_csv(genome_list_path, index=False)

    output_path = Path(output_path)
    for genome, (organism_id, genome_data) in cur_organism.items():
        genome_path = output_path / genome
        genome_path.mkdir(parents=True, exist_ok=True)

        download_file(
            f"https://imodulondb.org/organisms/{organism_id}/annotation/gene_files/gene_info.csv",
            genome_path / "gene_info.csv"
        )

        for imodulon_id, imodulon_name in genome_data:
            imodulon_path = genome_path / imodulon_id
            imodulon_path.mkdir(parents=True, exist_ok=True)

            download_file(
                f"https://imodulondb.org/organisms/{organism_id}/{imodulon_id}/data_files/gene_presence_list.csv",
                imodulon_path / "gene_presence_list.csv"
            )
            download_file(
                f"https://imodulondb.org/organisms/{organism_id}/{imodulon_id}/data_files/iM_table.csv",
                imodulon_path / "iM_table.csv"
            )
        

def run(args):
    organism_data(
        args.name,
        args.samples,
        args.genome_list,
        args.output,
    )
