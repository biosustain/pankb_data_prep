import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import json
from pathlib import Path
import argparse
import logging
import gzip

def initialize_parser(parser):
    parser.description = "Process gene-locustag relation info."
    parser.add_argument(
        "name",
        type=str,
        help="Family or analysis name to output data for.",
    )
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
        help="Directory containing per-gene folders with .fna and .faa files.",
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
        help="Output jsonl.gz file.",
    )

def gene_info(
    analysis_name, gp_locustag_path, all_locustag_path, fasta_dir, gtdb_meta_path, output_path
):
    fasta_dir = Path(fasta_dir)
    df_gene_presence_locustag = pd.read_csv(
        gp_locustag_path, index_col="Gene", low_memory=False
    )
    all_locustag_df = pd.read_csv(all_locustag_path, index_col=0, low_memory=False)
    all_locustag_df.rename(
        columns={
            "Locus_Tag": "locus_tag",
            "Genome_ID": "genome_id",
            "Gene_ID": "gene",
            "Prokka_Annotation": "protein",
            "Start_Position": "start_position",
            "End_Position": "end_position",
            "Nucleotide_Seq": "nucleotide_seq",
            "Amino_Acid_Seq": "aminoacid_seq",
        },
        inplace=True,
    )
    all_locustag_df.drop(["Nucleotide_Len", "less_than_2std_less_than_mean"], inplace=True, axis=1)
    df_gtdb_meta = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0)

    species = str(df_gtdb_meta.loc[df_gtdb_meta.index[0], "Organism"]).replace(
        "s__", ""
    )

    with gzip.open(output_path, "wt") as f:
    # with open(output_path, "w") as f:
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
                            logging.warning(
                                f"Locustag not in all_locustag_df: {genome_locustag}"
                            )
                else:
                    genome_locustag = f"{genome_id}@{locus_tag}"
                    if genome_locustag in all_locustag_df.index:
                        gene_locustag.append(genome_locustag)
                    else:
                        logging.warning(
                            f"Locustag not in all_locustag_df: {genome_locustag}"
                        )

            s_df = all_locustag_df.loc[gene_locustag, :].copy()
            s_df["nucleotide_seq"] = ""
            s_df["aminoacid_seq"] = ""
            s_df["species"] = species
            s_df["pangenome_analysis"] = analysis_name
            gene_locustag_only = [s.split("@", 1)[1] for s in gene_locustag]

            for k, ext in [("nucleotide_seq", "fna"), ("aminoacid_seq", "faa")]:
                for basename in ["pan_genes", "others"]:
                    gene_fasta_path = fasta_dir / gene_id / f"{basename}.{ext}"
                    if gene_fasta_path.is_file():
                        for record in SeqIO.parse(gene_fasta_path, "fasta"):
                            if record.id in gene_locustag_only:
                                lt = gene_locustag[gene_locustag_only.index(record.id)]
                                s = str(record.seq)
                                s_df.loc[lt, k] = s

            for record in s_df.to_dict(orient="records"):
                json.dump(
                    record,
                    f,
                    separators=(",", ":"),
                    ensure_ascii=False,
                    indent=None,
                )
                f.write("\n")


def run(args):
    gene_info(
        args.name,
        args.gp_locustag,
        args.all_locustag,
        args.fasta_dir,
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
