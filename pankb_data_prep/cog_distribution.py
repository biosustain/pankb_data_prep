import argparse
import json
import pandas as pd
from .utilities import COG_TABLE


def initialize_parser(parser):
    parser.description = "Generate COG distribution file."
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
        "--output_json",
        type=str,
        required=True,
        help="Output in json format.",
    )


def generate_cog_distribution(
    summary_v2_path, eggnog_summary_path, cog_distribution_path
):
    summary = pd.read_csv(summary_v2_path, index_col=0, low_memory=False)
    annotation = pd.read_csv(eggnog_summary_path, index_col=0, low_memory=False)
    annotation["COG_Categories"] = (
        "[" + annotation["COG_category"] + "]" + annotation["COG_category_name"]
    )
    cog_pan_class = pd.merge(
        summary["pangenome_class_2"],
        annotation[["COG_Categories", "COG_category_name", "COG_category"]],
        on=["Gene"],
    )

    # initialize the 3 classes list
    core = []
    accessory = []
    rare = []

    cog = COG_TABLE.Categories.values.tolist()

    for i in cog:
        core.append(
            len(
                cog_pan_class.loc[
                    (cog_pan_class["COG_category"].str.contains(i))
                    & (cog_pan_class["pangenome_class_2"]).str.contains("Core"),
                    :,
                ]
            )
        )
        accessory.append(
            len(
                cog_pan_class.loc[
                    (cog_pan_class["COG_category"].str.contains(i))
                    & (cog_pan_class["pangenome_class_2"]).str.contains("Accessory"),
                    :,
                ]
            )
        )
        rare.append(
            len(
                cog_pan_class.loc[
                    (cog_pan_class["COG_category"].str.contains(i))
                    & (cog_pan_class["pangenome_class_2"]).str.contains("Rare"),
                    :,
                ]
            )
        )

    geneclass_df = pd.DataFrame(
        {
            "Core": core,
            "Accessory": accessory,
            "Rare": rare,
            "Category": "[" + COG_TABLE.iloc[:, 1] + "]" + " " + COG_TABLE.iloc[:, 0],
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

    with open(cog_distribution_path, "w") as file:
        json.dump(class_json, file)


def run(args):
    generate_cog_distribution(
        args.summary,
        args.eggnog_summary,
        args.output_json,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
