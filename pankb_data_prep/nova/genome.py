import json
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import gzip

def initialize_parser(parser):
    parser.description = "Process data required for genome pages."
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
        "--gtdb_meta",
        type=str,
        required=True,
        help="GTDB meta csv file.",
    )
    parser.add_argument(
        "--mash_list",
        type=str,
        required=True,
        help="Mash list file.",
    )
    parser.add_argument(
        "--imodulon_dir",
        type=str,
        required=False,
        default=None,
        help="Directory contain iModulonDB info."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )
    parser.add_argument(
        "--iso_output",
        "-s",
        type=str,
        required=True,
        help="Output file or directory.",
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
        imodulons = [x.split(",", 1) for x in info_txt[1:]]
        s[genome] = (organism_id, imodulons)
    return s

def genome_info(
    analysis_name,
    species_summary_path,
    isosource_path,
    species_info_path,
    gp_binary_path,
    summary_v2_path,
    gtdb_meta_path,
    mash_list_path,
    imodulon_dir_path,
    genome_output_path,
    iso_output_path,
):
    genome_summary = pd.read_csv(species_summary_path, index_col=0, low_memory=False)
    isolation_src = pd.read_csv(isosource_path, index_col=0, low_memory=False)
    species_info = pd.read_csv(species_info_path, index_col=0, low_memory=False)
    phylo_group = pd.read_csv(mash_list_path, index_col=0, low_memory=False)
    phylo_group.rename(columns={"cluster": "phylo_group"}, inplace=True)

    if not "full_name" in species_info.columns:
        species_info["full_name"] = (species_info["genus"] + " " + species_info["species"] + " " + species_info["strain"]).str.strip()

    apm_binary = pd.read_csv(gp_binary_path, index_col=0, low_memory=False)
    summary = pd.read_csv(summary_v2_path, index_col=0, low_memory=False)
    genome_id_list = apm_binary.columns

    species_selection = set(genome_id_list)
    genome_info = pd.concat(
        [
            isolation_src.loc[list(species_selection & set(isolation_src.index)), :],
            genome_summary.loc[
                list(species_selection & set(genome_summary.index)),
                ["source", "gc_content", "genome_len"],
            ],
            species_info.loc[list(species_selection & set(species_info.index)), "full_name"],
            phylo_group.loc[list(species_selection & set(phylo_group.index)), "phylo_group"]
        ],
        axis=1,
    )

    imodulon_structure = get_imodulon_structure(imodulon_dir_path)

    genome_info.rename(columns={"Country": "country", "full_name": "strain"}, inplace=True)
    genome_info.drop(["biosample_accession", "source"], axis=1, inplace=True)

    df_gtdb_meta = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0)
    gtdb_meta_info = None
    for ind, info in df_gtdb_meta.iterrows():
        if info["Organism"] != "s__":
            gtdb_meta_info = info
            break
    species = str(gtdb_meta_info["Organism"]).replace(
        "s__", ""
    )
    family = str(gtdb_meta_info["Family"]).replace("f__", "")

    # Loop for all genome in the species
    # Get COG distribution in one genome
    with open(genome_output_path, "w") as f_genome:
        with open(iso_output_path, "w") as f_iso:
            for genome_id in genome_id_list:
                presence_gene_list = list(
                    (apm_binary.loc[apm_binary[genome_id] == 1, genome_id]).index
                )
                gene_class_distribution = [
                    int((summary.loc[presence_gene_list, "pangenome_class_2"] == pclass).sum())
                    for pclass in ["Core", "Accessory", "Rare"]
                ]
                genome_info_df = genome_info.loc[genome_id, :].copy()
                iso_info_df = genome_info.loc[genome_id, ["country", "geo_loc_name", "isolation_source"]].copy()
                iso_info_df["genome_id"] = genome_id
                genome_info_df.drop(["country", "geo_loc_name", "isolation_source"], inplace=True)
                genome_info_df["genome_id"] = genome_id
                genome_info_df["pangenome_analysis"] = analysis_name
                genome_info_df["species"] = species
                genome_info_df["gene_class_distribution"] = gene_class_distribution
                gi_index = genome_info_df.index.tolist()
                genome_info_df = genome_info_df.reindex(gi_index[-4:] + gi_index[:-4])

                if genome_id in imodulon_structure:
                    organism_id, imodulons = imodulon_structure[genome_id]
                    genome_info_df["imodulon_organism"] = organism_id
                    genome_info_df["imodulon_datasets"] = imodulons

                genome_record = genome_info_df.to_dict()
                json.dump(
                    genome_record,
                    f_genome,
                    separators=(",", ":"),
                    ensure_ascii=False,
                    indent=None,
                )
                f_genome.write("\n")

                iso_record = iso_info_df.to_dict()
                json.dump(
                    iso_record,
                    f_iso,
                    separators=(",", ":"),
                    ensure_ascii=False,
                    indent=None,
                )
                f_iso.write("\n")


def run(args):
    genome_info(
        args.name,
        args.species_summary,
        args.isosource,
        args.species_info,
        args.gp_binary,
        args.summary,
        args.gtdb_meta,
        args.mash_list,
        args.imodulon_dir,
        args.output,
        args.iso_output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
