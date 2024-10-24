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
        "--locus_tag_mapping",
        type=str,
        required=True,
        help="Locus tag mapping csv file.",
    )
    parser.add_argument(
        "--imodulon_tag_mapping",
        type=str,
        required=False,
        default=None,
        help="Locus tag mapping csv file.",
    )
    parser.add_argument(
        "--imodulon_dir",
        type=str,
        required=False,
        default=None,
        help="iModulonDB info dir.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output jsonl.gz file.",
    )

def get_imodulon_structure(imodulon_dir_path):
    s = {}
    if imodulon_dir_path is None:
        return s
    imodulon_dir_path = Path(imodulon_dir_path)
    if not imodulon_dir_path.is_dir():
        return s
    imodulon_genomes = [x.name for x in imodulon_dir_path.iterdir() if x.is_dir()]
    for genome in imodulon_genomes:
        genome_path = imodulon_dir_path / genome
        with open(genome_path / "info.txt", "r") as f:
            info_txt = [l_stripped for l in f.readlines() if (l_stripped := l.strip())]
        organism_id = info_txt[0]
        imodulons = []
        for x in info_txt[1:]:
            imodulon_dataset_id, imodulon_dataset_name = x.split(",", 1)
            df_iM_gene_presence = pd.read_csv(genome_path / imodulon_dataset_id / "gene_presence_list.csv", index_col="Gene")
            df_iM_table = pd.read_csv(genome_path / imodulon_dataset_id / "iM_table.csv", index_col="k")            
            imodulons.append((imodulon_dataset_id, imodulon_dataset_name, df_iM_gene_presence, df_iM_table))
        s[genome] = (organism_id, imodulons)
    return s

def gene_info(
    analysis_name, gp_locustag_path, all_locustag_path, fasta_dir, gtdb_meta_path, locus_tag_mapping_path, imodulon_tag_mapping_path, imodulon_dir_path, output_path
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

    df_locus_tag_mapping = pd.read_csv(locus_tag_mapping_path, index_col=["prokka_locus_tag","genome"])
    if imodulon_tag_mapping_path is None:
        df_imodulon_tag_mapping = None
    else:
        df_imodulon_tag_mapping = pd.read_csv(imodulon_tag_mapping_path, index_col=["prokka_locus_tag","genome"])

    species = str(df_gtdb_meta.loc[df_gtdb_meta.index[0], "Organism"]).replace(
        "s__", ""
    )

    imodulon_structure = get_imodulon_structure(imodulon_dir_path)

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

            for ind, row in s_df.items():
                try:
                    lt_map = df_locus_tag_mapping.loc[(row["locus_tag"], row["genome_id"]), :]
                except:
                    continue
                row["original_locus_tag"] = lt_map["original_locus_tag"]
                row["original_gene"] = lt_map["original_gene"]
                row["original_exact_match"] = lt_map["exact_match"]
            
            if not df_imodulon_tag_mapping is None:
                for ind, row in s_df.items():
                    if not row["genome_id"] in imodulon_structure:
                        continue
                    try:
                        lt_map = df_imodulon_tag_mapping.loc[(row["locus_tag"], row["genome_id"]), :]
                    except:
                        continue
                    iM_locus_tag = lt_map["original_locus_tag"]

                    organism_id, imodulon_datasets = imodulon_structure[row["genome_id"]]
                    d = []
                    for imodulon_dataset_id, imodulon_dataset_name, df_iM_gene_presence, df_iM_table in imodulon_datasets:
                        if not iM_locus_tag in df_iM_gene_presence.index:
                            continue
                        iM_k = df_iM_gene_presence.loc[iM_locus_tag, "iModulon"]
                        iM_name = df_iM_table.loc[iM_k, "name"]
                        iM_data = {
                            "locus_tag": iM_locus_tag,
                            "organism": organism_id,
                            "dataset": imodulon_dataset_id,
                            "dataset_name": imodulon_dataset_name,
                            "k": iM_k,
                            "name": iM_name,
                            "exact_match": lt_map["exact_match"],
                        }
                        d.append(iM_data)
                    row["imodulon_data"] = d

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
        args.locus_tag_mapping,
        args.imodulon_tag_mapping,
        args.imodulon_dir,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
