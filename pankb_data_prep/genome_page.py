import json
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from .utilities import COG_TABLE


def initialize_parser(parser):
    parser.description = "Process data required for genome pages."
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
        "--eggnog_summary",
        type=str,
        required=True,
        help="Eggnog summary file.",
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
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )


def generate_genome_page(
    species_summary_path,
    isosource_path,
    species_info_path,
    gp_binary_path,
    summary_v2_path,
    eggnog_summary_path,
    genome_page_dir,
):
    genome_page_dir = Path(genome_page_dir)
    genome_summary = pd.read_csv(species_summary_path, index_col=0, low_memory=False)
    isolation_src = pd.read_csv(isosource_path, index_col=0, low_memory=False)
    species_info = pd.read_csv(species_info_path, index_col=0, low_memory=False)
    # species_info["full_name"] = (
    #     species_info.genus + " " + species_info.species + " " + species_info.strain
    # )
    # species_info_selected = species_info.loc[list(isolation_src.index), :]
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

    apm_binary = pd.read_csv(gp_binary_path, index_col=0, low_memory=False)
    summary = pd.read_csv(summary_v2_path, index_col=0, low_memory=False)
    annotation = pd.read_csv(eggnog_summary_path, index_col=0, low_memory=False)

    annotation["COG_Categories"] = (
        "[" + annotation["COG_category_name"] + "]" + annotation["COG_category"]
    )
    cog_pan_class = pd.merge(
        summary["pangenome_class_2"],
        annotation[["COG_Categories", "COG_category_name", "COG_category"]],
        on=["Gene"],
    )
    genome_id_list = apm_binary.columns

    # Loop for all genome in the species
    # Get COG distribution in one genome
    for genome_id in genome_id_list:

        # initialize the 3 classes list
        core = []
        accessory = []
        rare = []

        cog = COG_TABLE.Categories.values.tolist()
        presence_gene_list = list(
            (apm_binary.loc[apm_binary[genome_id] == 1, genome_id]).index
        )
        cog_presence = cog_pan_class.loc[presence_gene_list, :]

        for i in cog:
            core.append(
                len(
                    cog_presence.loc[
                        (cog_presence["COG_category"].str.contains(i))
                        & (cog_presence["pangenome_class_2"]).str.contains("Core"),
                        :,
                    ]
                )
            )
            accessory.append(
                len(
                    cog_presence.loc[
                        (cog_presence["COG_category"].str.contains(i))
                        & (cog_presence["pangenome_class_2"]).str.contains("Accessory"),
                        :,
                    ]
                )
            )
            rare.append(
                len(
                    cog_presence.loc[
                        (cog_presence["COG_category"].str.contains(i))
                        & (cog_presence["pangenome_class_2"]).str.contains("Rare"),
                        :,
                    ]
                )
            )

        geneclass_df = pd.DataFrame(
            {
                "Core": core,
                "Accessory": accessory,
                "Rare": rare,
                "Category": "["
                + COG_TABLE.iloc[:, 0]
                + "]"
                + " "
                + COG_TABLE.iloc[:, 1],
            }
        )

        geneclass_df["All"] = geneclass_df.iloc[:, :3].sum(axis=1)
        geneclass_df_sorted = geneclass_df.sort_values(by=["All"], ascending=False)

        class_json = {
            "categories": geneclass_df_sorted["Category"].values.tolist(),
            "Core": geneclass_df_sorted["Core"].values.tolist(),
            "Accessory": geneclass_df_sorted["Accessory"].values.tolist(),
            "Rare": geneclass_df_sorted["Rare"].values.tolist(),
        }

        genome_id_folder = genome_page_dir / genome_id
        genome_id_folder.mkdir(parents=True, exist_ok=True)

        cog_distr_file = genome_id_folder / "COG_distribution.json"

        with open(cog_distr_file, "w") as file:
            json.dump(class_json, file)

        # Genome info
        genome_info_file = genome_id_folder / "genome_info.json"
        gene_class_distribution = [
            sum(cog_presence["pangenome_class_2"] == "Core"),
            sum(cog_presence["pangenome_class_2"] == "Accessory"),
            sum(cog_presence["pangenome_class_2"] == "Rare"),
        ]
        genome_info_df = pd.DataFrame(
            [genome_info.loc[genome_id, :]], columns=genome_info.columns
        )
        genome_info_df["Gene_class_distribution"] = [gene_class_distribution]
        genome_info_df.to_json(path_or_buf=genome_info_file, orient="index")


def run(args):
    generate_genome_page(
        args.species_summary,
        args.isosource,
        args.species_info,
        args.gp_binary,
        args.summary,
        args.eggnog_summary,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
