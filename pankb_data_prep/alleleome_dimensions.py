import pandas as pd
import argparse


def initialize_parser(parser):
    parser.description = "Generate full summary table."
    parser.add_argument(
        "--sel_genes",
        type=str,
        required=True,
        help="alleleome sel_genes csv file.",
    )
    parser.add_argument(
        "--filt_norm",
        type=str,
        required=True,
        help="Alleleome filt_norm file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )


def alleleome_dims(sel_genes_path, filt_norm_path, dims_path):
    sel_genes = pd.read_csv(sel_genes_path, low_memory=False, index_col=0)

    n_genes = int(sel_genes["Pan"].values.sum())

    df_filt_norm = pd.read_csv(filt_norm_path, header=0, index_col=0, low_memory=False)

    n_muts = int((df_filt_norm["Sequence_type"] == "Variant").values.sum())

    df = pd.DataFrame({"N_genes": [n_genes], "N_muts": [n_muts]}, index=["Unknown"])
    df.index.name = "Project"

    df.to_csv(dims_path)


def run(args):
    alleleome_dims(args.sel_genes, args.filt_norm, args.output)


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
