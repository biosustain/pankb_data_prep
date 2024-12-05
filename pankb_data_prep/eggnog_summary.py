import logging
from pathlib import Path
import math
import pandas as pd
from Bio import SeqIO
import argparse
from .utilities import COG_DICT


def initialize_parser(parser):
    parser.description = "Parse Eggnog annotations to a summary file."
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
        "--eggnog_table",
        type=str,
        required=True,
        help="Eggnog table file (.emapper.annotations).",
    )
    parser.add_argument(
        "--reference",
        type=str,
        required=True,
        help="Pangenome reference fasta file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )

def generate_eggnog_summary(
    gp_binary_path, summary_path, eggnog_table_path, reference_path, eggnog_summary_path
):
    df_roary_binary = pd.read_csv(gp_binary_path, index_col=0, low_memory=False)
    df_pangene_summary = pd.read_csv(summary_path, index_col=0)
    # eggnog_table_path = '/home/binhuan/data_pankb/bgcflow/data/interim/eggnog_roary/' + species + '/' + species + '.emapper.annotations'

    df_eggnog = pd.read_csv(
        eggnog_table_path, sep="\t", header=4, index_col="#query"
    ).iloc[:-3, :]
    df_eggnog.index.name = "locus_tag"

    # Sort the matrix by the sum of strains presence
    df = pd.DataFrame()
    df["data"] = df_roary_binary.sum(axis=1).sort_values(ascending=True)

    # Group the data by value and count the number of occurrences
    df = df.groupby("data").size().reset_index(name="frequency")
    x15 = math.floor(df["data"].quantile(0.15))
    x99 = math.floor(df["data"].quantile(0.99).round()) - 1

    df_accesory = df_roary_binary[
        (df_roary_binary.sum(1) >= x15) & (df_roary_binary.sum(1) < x99)
    ]

    # fasta_file = '/home/binhuan/data_pankb/bgcflow/data/interim/roary/' + species + '/pan_genome_reference.fa'
    recs = SeqIO.parse(reference_path, format="fasta")
    for seq in recs:
        desc = seq.description
        locus_tag = desc.split(" ")[0]
        group_id = desc.split(" ")[1]
        if locus_tag in df_eggnog.index:
            df_eggnog.loc[locus_tag, "group_id"] = group_id

    for pan_gene_id in df_pangene_summary.index:
        if pan_gene_id in df_eggnog.group_id.tolist():
            locus_tag = df_eggnog[df_eggnog.group_id == pan_gene_id].index.tolist()[0]
            for col in df_eggnog.columns:
                df_pangene_summary.loc[pan_gene_id, col] = df_eggnog.loc[locus_tag, col]

            cog_cat = df_eggnog.loc[locus_tag, "COG_category"]
            if len(cog_cat) > 1:
                cog_name = " | ".join([COG_DICT[cog_letter] for cog_letter in cog_cat])
                df_pangene_summary.loc[pan_gene_id, "COG_category_name"] = cog_name
                df_pangene_summary.loc[pan_gene_id, "uniq_COG_category_name"] = (
                    "multi_COGs_" + str(len(cog_cat))
                )
            else:
                df_pangene_summary.loc[pan_gene_id, "COG_category_name"] = COG_DICT[
                    cog_cat
                ]
                df_pangene_summary.loc[pan_gene_id, "uniq_COG_category_name"] = (
                    COG_DICT[cog_cat]
                )
        else:
            for col in df_eggnog.columns:
                df_pangene_summary.loc[pan_gene_id, col] = "-"
            df_pangene_summary.loc[pan_gene_id, "COG_category_name"] = (
                "Not found in COG"
            )
            df_pangene_summary.loc[pan_gene_id, "uniq_COG_category_name"] = (
                "Not found in COG"
            )

    df_pangene_summary.to_csv(eggnog_summary_path)
    # df_pangene_summary.to_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_eggnog_summary.csv')


def run(args):
    generate_eggnog_summary(
        args.gp_binary,
        args.summary,
        args.eggnog_table,
        args.reference,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
