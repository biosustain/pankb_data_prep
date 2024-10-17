import argparse
import pandas as pd
import json
from ..utilities import calculate_lambda

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
        "--filt_norm",
        type=str,
        required=True,
        help="Alleleome filt_norm csv file.",
    )
    parser.add_argument(
        "--gtdb_meta",
        type=str,
        required=True,
        help="GTDB meta csv file.",
    )
    parser.add_argument(
        "--sel_genes",
        type=str,
        required=True,
        help="alleleome sel_genes csv file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )

def openness_check(mean_lambda):
    if mean_lambda < 0.3:
        return 'Closed'
    else:
        return 'Intermediate Open'

def organism_info(
    analysis_name,
    gp_binary_path,
    summary_v2_path,
    gtdb_meta_path,
    filt_norm_path,
    sel_genes_path,
    output_path,
):
    df_gp_binary = pd.read_csv(gp_binary_path, index_col="Gene", low_memory=False)
    df_pangene_summary = pd.read_csv(summary_v2_path, low_memory=False, index_col=0)
    df_gtdb_meta = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0)
    df_filt_norm = pd.read_csv(filt_norm_path, header=0, index_col=False, usecols=['Gene', 'Sequence_type'])
    df_filt_norm = df_filt_norm.loc[df_filt_norm["Sequence_type"] == "Variant", "Gene"]
    sel_genes = pd.read_csv(sel_genes_path, low_memory=False, index_col=0)

    genomes = df_gp_binary.columns.tolist()
    n_genomes = len(genomes)
    n_genes = int(df_pangene_summary.shape[0])
    n_alleleomes = int(sel_genes["Pan"].sum())
    n_muts = int(df_filt_norm.shape[0])

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

    json_data = {
        "family": family,
        "species": species,
        "pangenome_analysis": analysis_name,
        "genomes_num": n_genomes,
        "gene_class_distribution": [core_len, accessory_len, rare_len],
        "openness": openness,
        "genes_num": n_genes,
        "alleleomes_num": n_alleleomes,
        "mutations_num": n_muts,
    }
    with open(output_path, "w") as f:
        json.dump(json_data, f)


def run(args):
    organism_info(
        args.name,
        args.gp_binary,
        args.summary,
        args.gtdb_meta,
        args.filt_norm,
        args.sel_genes,
        args.output,
    )
