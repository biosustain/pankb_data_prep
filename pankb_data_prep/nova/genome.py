import json
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import gzip

def initialize_parser(parser):
    parser.description = "Process data required for genome pages."
    parser.add_argument(
        "name",
        type=str,
        help="Family or analysis name to output data for.",
    )
    parser.add_argument(
        "--gp_binary",
        type=str,
        required=True,
        help="Gene presence binary csv file.",
    )
    parser.add_argument(
        "--summary",
        type=str,
        required=True,
        help="Pangene summary csv file.",
    )
    parser.add_argument(
        "--species_summary",
        type=str,
        required=True,
        help="Species summary csv file.",
    )
    parser.add_argument(
        "--isosource",
        type=str,
        required=True,
        help="Isolation source file.",
    )
    parser.add_argument(
        "--species_info",
        type=str,
        required=True,
        help="Species info df_ncbi_meta csv file.",
    )
    parser.add_argument(
        "--gtdb_meta",
        type=str,
        required=True,
        help="GTDB meta csv file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )


def genome_info(
    analysis_name,
    species_summary_path,
    isosource_path,
    species_info_path,
    gp_binary_path,
    summary_v2_path,
    gtdb_meta_path,
    output_path,
):
    genome_summary = pd.read_csv(species_summary_path, index_col=0, low_memory=False)
    isolation_src = pd.read_csv(isosource_path, index_col=0, low_memory=False)
    species_info = pd.read_csv(species_info_path, index_col=0, low_memory=False)
    species_selection = list(isolation_src.index)
    genome_info = pd.concat(
        [
            isolation_src.loc[species_selection, :],
                genome_summary.loc[
                    species_selection,
                    ["source", "gc_content", "genome_len"],
                ],
                species_info.loc[species_selection, "full_name"]
        ],
        axis=1,
    )

    genome_info.rename(columns={"Country": "country", "full_name": "strain"}, inplace=True)
    genome_info.drop(["biosample_accession", "source"], axis=1, inplace=True)

    apm_binary = pd.read_csv(gp_binary_path, index_col=0, low_memory=False)
    summary = pd.read_csv(summary_v2_path, index_col=0, low_memory=False)
    genome_id_list = apm_binary.columns

    df_gtdb_meta = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0)
    species = str(df_gtdb_meta.loc[df_gtdb_meta.index[0], "Organism"]).replace(
        "s__", ""
    )

    # Loop for all genome in the species
    # Get COG distribution in one genome
    with open(output_path, "w") as f:
        for genome_id in genome_id_list:
            presence_gene_list = list(
                (apm_binary.loc[apm_binary[genome_id] == 1, genome_id]).index
            )
            gene_class_distribution = [
                int((summary.loc[presence_gene_list, "pangenome_class_2"] == pclass).sum())
                for pclass in ["Core", "Accessory", "Rare"]
            ]
            genome_info_df = genome_info.loc[genome_id, :].copy()
            genome_info_df["genome_id"] = genome_id
            genome_info_df["pangenome_analysis"] = analysis_name
            genome_info_df["species"] = species
            genome_info_df["gene_class_distribution"] = gene_class_distribution
            gi_index = genome_info_df.index.tolist()
            genome_info_df = genome_info_df.reindex(gi_index[-4:] + gi_index[:-4])
            record = genome_info_df.to_dict()
            json.dump(
                record,
                f,
                separators=(",", ":"),
                ensure_ascii=False,
                indent=None,
            )
            f.write("\n")


def run(args):
    genome_info(
        args.name,
        args.species_summary,
        args.isosource,
        args.species_info,
        args.gp_binary,
        args.summary,
        args.gtdb_meta,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
