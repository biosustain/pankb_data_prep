import argparse
import pandas as pd
import os
import json
import requests
from io import StringIO
from ratelimiter import RateLimiter


def initialize_parser(parser):
    parser.description = "Process pangene data."
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
        "--eggnog_table",
        type=str,
        required=True,
        help="Eggnog table file.",
    )
    parser.add_argument(
        "--summary",
        type=str,
        required=True,
        help="Pangene summary v2 csv file.",
    )
    parser.add_argument(
        "--gtdb_meta",
        type=str,
        required=True,
        help="GTDB meta csv file.",
    )
    parser.add_argument(
        "--pangene_output",
        "-o",
        type=str,
        required=True,
        help="Output in json format.",
    )
    parser.add_argument(
        "--pathway_output",
        "-p",
        type=str,
        required=True,
        help="Output in json format.",
    )


COG_DICT = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport, and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
    "-": "Not found in COG",
}


def expand_cog_data(row):
    cog_cat = row["COG_category"]
    if pd.isna(cog_cat):
        cog_cat = "-"
    cog_name = " | ".join([COG_DICT[cog_letter] for cog_letter in cog_cat])
    return [
        cog_cat,
        cog_name,
        (f"multi_COGs_{len(cog_cat):d}" if len(cog_cat) > 0 else cog_cat),
    ]


@RateLimiter(max_calls=2, period=1)
def get_kegg_info(query):
    url = f"https://rest.kegg.jp/list/{query}"
    print(url)
    response = requests.get(url)
    df = pd.read_csv(
        StringIO(response.text),
        names=["id", "name"],
        index_col=0,
        sep="\t",
        header=None,
        skipinitialspace=True,
    )
    return df


def pangene_info(
    analysis_name,
    eggnog_table_path,
    summary_v2_path,
    gp_binary_path,
    gtdb_meta_path,
    pangene_output_path,
    pathway_output_path,
):
    df_summary = pd.read_csv(summary_v2_path, low_memory=False, index_col=0)
    gp_binary = pd.read_csv(gp_binary_path, low_memory=False, index_col=0)
    gp_binary["Occurency"] = gp_binary.iloc[:, 1:].sum(axis=1)
    df_gtdb_meta = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0)

    df_eggnog = pd.read_csv(
        eggnog_table_path,
        sep="\t",
        header=4,
        index_col="#query",  # na_values="-",
        dtype={"eggNOG_OGs": str},
    ).iloc[:-3, :]
    df_eggnog.index.name = "locus_tag"

    df = pd.merge(df_summary, gp_binary[["Occurency"]], on="Gene", how="left")
    df.reset_index(inplace=True)
    df = pd.merge(df, df_eggnog, on="locus_tag", how="left")
    df.set_index("Gene", inplace=True)

    species = str(df_gtdb_meta.loc[df_gtdb_meta.index[0], "Organism"]).replace(
        "s__", ""
    )
    family = str(df_gtdb_meta.loc[df_gtdb_meta.index[0], "Family"]).replace("f__", "")

    df[["COG_category", "COG_category_name", "uniq_COG_category_name"]] = df.apply(
        expand_cog_data, axis=1, result_type="expand"
    )
    df.fillna("-", inplace=True)

    key_mapping = {
        v: v.lower()
        for v in [
            "PFAMs",
            "EC",
            "eggNOG_OGs",
            "GOs",
            "KEGG_ko",
            "KEGG_Module",
            "KEGG_Reaction",
            "KEGG_TC",
            "BRITE",
            "CAZy",
        ]
    }

    gene_data = []
    pathway_data = {}
    for gene, row in df.iterrows():
        pathways = row["KEGG_Pathway"].split(",") if row["KEGG_Pathway"] != "-" else []
        pathways = [
            p
            for p in pathways
            if (p.startswith("map") or (not p.replace("ko", "map") in pathways))
        ]  # Prefer the map links over ko
        entry = {
            "gene": gene,
            "pangenome_analysis": analysis_name,
            "species": species,
            "family": family,
            "cog_category": row["COG_category"],
            "cog_name": row["COG_category_name"],
            "description": row["Description"],
            "protein": row["Annotation"],
            "frequency": row["Occurency"],
            "pangenomic_class": row["pangenome_class_2"],
        }
        if pathways:
            entry["kegg_pathway"] = pathways
        for k, v in key_mapping.items():
            if row[k] != "-":
                entry[v] = row[k].split(",")
        gene_data.append(entry)

        for pathway in pathways:
            pathway_data[pathway] = pathway_data.get(pathway, []) + [
                {"pangenome_analysis": analysis_name, "gene": gene}
            ]

    with open(pangene_output_path, "w") as f:
        json.dump(gene_data, f)

    kegg_pathway_names = get_kegg_info("pathway")
    q = []
    pwdata_len = len(pathway_data) - 1
    for i, pathway in enumerate(pathway_data.keys()):
        if not pathway in kegg_pathway_names.index:
            q.append(pathway)
        if len(q) == 10 or i == pwdata_len:
            new_pathway_names = get_kegg_info("+".join(q))
            kegg_pathway_names = pd.concat([kegg_pathway_names, new_pathway_names])
            q = []
    kegg_pathway_names = kegg_pathway_names["name"]

    pathway_data = [
        {
            "pathway_id": k,
            "pathway_name": kegg_pathway_names.get(k, "Unknown"),
            "genes": v,
        }
        for k, v in pathway_data.items()
    ]

    with open(pathway_output_path, "w") as f:
        json.dump(pathway_data, f)


def run(args):
    pangene_info(
        args.name,
        args.eggnog_table,
        args.summary,
        args.gp_binary,
        args.gtdb_meta,
        args.pangene_output,
        args.pathway_output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
