import argparse
import pandas as pd
import json
from .utilities import calculate_lambda

def initialize_parser(parser):
    parser.description = "Generate species summary."
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
    parser.add_argument(
        "--output_json",
        type=str,
        required=True,
        help="Output in json format.",
    )

def openness_check(mean_lambda):
    if mean_lambda < 0.3:
        return 'Closed'
    else:
        return 'Intermediate Open'

def species_pangenome_summary(
    analysis_name,
    gp_binary_path,
    summary_v2_path,
    gtdb_meta_path,
    species_summary_csv_path,
    species_summary_json_path,
):
    df_gp_binary = pd.read_csv(gp_binary_path, index_col="Gene", low_memory=False)
    df_pangene_summary = pd.read_csv(summary_v2_path, low_memory=False, index_col=0)
    df_gtdb_meta = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0)

    genomes = df_gp_binary.columns.tolist()
    n_genomes = len(genomes)
    n_genes = int(df_pangene_summary.shape[0])

    # TODO
    # print(set(df_gtdb_meta.loc[genomes, "Family"]))
    # print(set(df_gtdb_meta.loc[genomes, "Species"]))
    # assert len(set(df_gtdb_meta.loc[genomes, "Family"])) == 1 # Make sure that all the genomes in the analysis are from the same family
    # assert len(set(df_gtdb_meta.loc[genomes, "Species"])) == 1 # Make sure that all the genomes in the analysis are from the same species

    species = str(df_gtdb_meta.loc[genomes[0], "Organism"]).replace("s__", "")
    family = str(df_gtdb_meta.loc[genomes[0], "Family"]).replace("f__", "")

    core_len = int(
        df_pangene_summary.loc[
            df_pangene_summary["pangenome_class_2"] == "Core", "pangenome_class_2"
        ].count()
    )
    rare_len = int(
        df_pangene_summary.loc[
            df_pangene_summary["pangenome_class_2"] == "Rare", "pangenome_class_2"
        ].count()
    )
    accessory_len = int(
        df_pangene_summary.loc[
            df_pangene_summary["pangenome_class_2"] == "Accessory", "pangenome_class_2"
        ].count()
    )

    #TODO Remove .T by changing calculate_lambda to sample different axis
    mean_lambda, std_lambda = calculate_lambda(df_gp_binary.T)
    openness = openness_check(mean_lambda)
    # openness = "Open"

    df = pd.DataFrame(
        {
            "Pangenome_analyses": [analysis_name],
            "Family": [family],
            "Species": [species],
            "N_of_genome": [n_genomes],
            "N_of_gene": [n_genes],
            "N_of_core": [core_len],
            "N_of_rare": [rare_len],
            "N_of_accessory": [accessory_len],
            "Openness": [openness],
        }
    )
    df.to_csv(species_summary_csv_path, index=False)

    json_data = {
        "Family": family,
        "Species": species,
        "Number_of_genome": n_genomes,
        "Gene_class": [core_len, accessory_len, rare_len],
        "Openness": openness,
    }
    with open(species_summary_json_path, "w") as f:
        json.dump(json_data, f)


def run(args):
    species_pangenome_summary(
        args.name,
        args.gp_binary,
        args.summary,
        args.gtdb_meta,
        args.output,
        args.output_json,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
