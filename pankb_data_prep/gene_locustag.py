import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import json
from pathlib import Path
import argparse
import logging


def initialize_parser(parser):
    parser.description = "Process gene-locustag relation info."
    parser.add_argument(
        "--gp_locustag",
        type=str,
        required=True,
        help="Gene presence locustag csv file.",
    )
    parser.add_argument(
        "--all_locustag",
        type=str,
        required=True,
        help="Locustag info df_all_locustag csv file.",
    )
    parser.add_argument(
        "--fasta_dir",
        type=str,
        required=True,
        help="Directory containing per-gene folders with .fna and .faa files."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )


def remove_slash(s):
    if "/" in str(s):
        return s.replace("/", "_")
    else:
        return s


def generate_locustag_data(gp_locustag_path, all_locustag_path, fasta_dir, gene_locustag_dir):
    gene_locustag_dir = Path(gene_locustag_dir)
    gene_locustag_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir = Path(fasta_dir)
    df_gene_presence_locustag = pd.read_csv(
        gp_locustag_path, index_col="Gene", low_memory=False
    )
    df_gene_presence_locustag.index = [
        remove_slash(i) for i in list(df_gene_presence_locustag.index)
    ]
    all_locustag_df = pd.read_csv(all_locustag_path, index_col=0, low_memory=False)

    for gene_id in df_gene_presence_locustag.index.tolist():
        gene_locustag = []
        locus_tag_series = df_gene_presence_locustag.loc[gene_id, :].dropna()
        locus_tag_series = locus_tag_series.str.split("\t")
        locus_tag_series = locus_tag_series.explode()
        for genome_id, locus_tag in locus_tag_series.items():
            if "\t" in locus_tag:
                locus_tag_list = locus_tag.split("\t")
                for locus_tag_id in locus_tag_list:
                    genome_locustag = f"{genome_id}@{locus_tag_id}"
                    if genome_locustag in all_locustag_df.index:
                        gene_locustag.append(genome_locustag)
                    else:
                        logging.warn(f"Locustag not in all_locustag_df: {gene_locustag}")
            else:
                genome_locustag = f"{genome_id}@{locus_tag}"
                if genome_locustag in all_locustag_df.index:
                    gene_locustag.append(genome_locustag)
                else:
                    logging.warn(f"Locustag not in all_locustag_df: {gene_locustag}")

        s_df = all_locustag_df.loc[gene_locustag, :].copy()
        s_df["Nucleotide_Seq"] = ""
        s_df["Amino_Acid_Seq"] = ""
        gene_locustag_only = [s.split("@", 1)[1] for s in gene_locustag]

        for k, ext in [("Nucleotide_Seq", "fna"), ("Amino_Acid_Seq", "faa")]:
            for basename in ["pan_genes", "others"]:
                gene_fasta_path = fasta_dir / gene_id / f"{basename}.{ext}"
                if gene_fasta_path.is_file():
                    for record in SeqIO.parse(gene_fasta_path, "fasta"):
                        if record.id in gene_locustag_only:
                            lt = gene_locustag[gene_locustag_only.index(record.id)]
                            s_df.loc[lt, k] = str(record.seq)


        # Write the JSON object to a file
        with open(gene_locustag_dir / f"{gene_id}.json", "w") as f:
            json.dump(
                s_df.to_dict(orient="records"),
                f,
                separators=(",", ":"),
                ensure_ascii=False,
            )


def run(args):
    generate_locustag_data(
        args.gp_locustag,
        args.all_locustag,
        args.fasta_dir,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
