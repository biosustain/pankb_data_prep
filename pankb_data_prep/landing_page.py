import pandas as pd
import json
import argparse


def initialize_parser(parser):
    parser.description = "Process data required for the landing page."
    parser.add_argument(
        "--species_summary",
        type=str,
        required=True,
        help="Species summary csv file.",
    )
    parser.add_argument(
        "--alleleome_dimensions",
        type=str,
        required=True,
        help="Alleleome dimensions csv file.",
    )
    parser.add_argument(
        "--output_pankb_dimension",
        type=str,
        required=True,
        help="Gene count json file.",
    )
    parser.add_argument(
        "--output_json",
        type=str,
        required=True,
        help="Output in json format.",
    )


def generate_landing_page(
    species_summary_path, alleleome_dimensions_path, pankb_dimension_path, species_genome_gene_path,
):
    df_species_summaries = pd.read_csv(species_summary_path, index_col=0)
    df_alleleome_dimensions = pd.read_csv(alleleome_dimensions_path, index_col=0)


    N_alleleome = int(df_alleleome_dimensions["N_genes"].sum())
    N_mutations = int(df_alleleome_dimensions["N_muts"].sum())

    gene_cluster_count = int(df_species_summaries["N_of_gene"].sum())
    genome_count = int(df_species_summaries["N_of_genome"].sum())
    species_count = int(df_species_summaries.shape[0])
    dimension = {
        "Mutations": N_mutations,
        "Genes": gene_cluster_count,
        "Alleleomes": N_alleleome,
        "Genomes": genome_count,
        "Species pangenomes": species_count,
    }
    with open(pankb_dimension_path, "w") as f:
        json.dump(dimension, f)

    family_genome_gene = {}
    for family in sorted(set(df_species_summaries["Family"].tolist())):
        family_genome_gene[family] = {}

        for species, row in df_species_summaries.loc[
            df_species_summaries["Family"] == family, :
        ].iterrows():
            family_genome_gene[family][species] = [
                int(row["N_of_genome"]),
                int(row["N_of_gene"]),
            ]

    with open(species_genome_gene_path, "w") as f:
        json.dump(family_genome_gene, f)


def run(args):
    generate_landing_page(
        args.species_summary,
        args.alleleome_dimensions,
        args.output_pankb_dimension,
        args.output_json,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
