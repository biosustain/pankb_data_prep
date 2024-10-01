import json
import pandas as pd
import numpy as np
import gzip
import logging
import argparse


def initialize_parser(parser):
    parser.description = "Process data for the heatmap."
    parser.add_argument(
        "--gp_locustag",
        type=str,
        required=True,
        help="Gene presence locustag csv file.",
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
        "--eggnog_summary",
        type=str,
        required=True,
        help="Eggnog summary file.",
    )
    parser.add_argument(
        "--mash_list",
        type=str,
        required=True,
        help="Mash list file.",
    )
    parser.add_argument(
        "--isosource",
        type=str,
        required=True,
        help="Isolation source file.",
    )
    # parser.add_argument(
    #     "--species_info",
    #     type=str,
    #     required=True,
    #     help="Species info df_ncbi_meta csv file.",
    # )
    parser.add_argument(
        "--heatmap_base",
        type=str,
        required=True,
        help="Basename for the heatmap output, class + '.json.gz' will be added.",
    )
    parser.add_argument(
        "--sourceinfo_base",
        type=str,
        required=True,
        help="Basename for the source info output, class + '.json' will be added.",
    )


def gzip_file(input_path, output_path):
    with open(input_path, "rb") as file_in:
        with gzip.open(output_path, "wb") as file_out:
            file_out.writelines(file_in)


def filter_cluster(l):
    top2 = pd.to_numeric(l).nlargest(2)
    if top2.iloc[0] >= 0.7 and top2.iloc[1] <= 0.3:
        return pd.to_numeric(l).idxmax()
    else:
        return False


def rows(source):
    source = source.replace("Missing", "undefined")

    # get the genomes' name (index label)
    genome = source.index

    genome_list = []
    for i in range(len(genome)):
        genome_list.append({"name": genome[i], "meta": [str(source.iloc[i, 4])]})
    return genome_list


def cols(summary, annotation, locustag_ordered):
    annotation["annotation_all"] = (
        "[" + annotation["COG_category"] + "]" + " " + annotation["COG_category_name"]
    )
    merged_df = pd.merge(
        annotation.loc[:, ["annotation_all", "Description"]],
        summary["pangenome_class_2"],
        on=["Gene"],
    )
    anno_list = merged_df.values.tolist()

    # get the number of genes and genomes and their name
    gene_num = len(anno_list)
    genome_num = len(locustag_ordered.columns)
    gene_name = merged_df.index
    genome_name = locustag_ordered.columns
    gene_list = []

    for i in range(gene_num):
        gene_list.append(
            {"name": gene_name[i], "meta": merged_df.iloc[i, :].values.tolist()}
        )
    return gene_list


def matrix(gp_binary):
    return gp_binary.values.tolist()


def source_info(source):
    source = source.replace("Missing", "undefined")

    # get the genomes' name (index label)
    genome = source.index

    genome_dict = {}
    for i in range(len(genome)):
        genome_dict[str(genome[i])] = [
            str(source.iloc[i, 1]),
            str(source.iloc[i, 0]),
            str(source.iloc[i, 3]),
        ]
    return genome_dict


def remove_slash(s):
    if "/" in str(s):
        return s.replace("/", "_")
    else:
        return s


def generate_heatmap(
    gp_binary_path,
    gp_locustag_path,
    summary_v2_path,
    eggnog_summary_path,
    mash_list_path,
    isosource_path,
    # species_info_path,
    heatmap_base_path,
    sourceinfo_base_path,
):
    # load data
    gp_binary = pd.read_csv(gp_binary_path, index_col=0, low_memory=False).T
    gp_binary.index.names = ["genome_id"]
    summary = pd.read_csv(summary_v2_path, index_col=0, low_memory=False)
    annotation = pd.read_csv(eggnog_summary_path, index_col=0, low_memory=False)
    gp_locustag = pd.read_csv(gp_locustag_path, index_col=0, low_memory=False)
    phylo_group = pd.read_csv(mash_list_path, index_col=0, low_memory=False)
    isolation_src = pd.read_csv(isosource_path, index_col=0, low_memory=False)
    # species_info = pd.read_csv(species_info_path, index_col=0, low_memory=False)

    # species_info["genome_name"] = (
    #     species_info["genus"]
    #     + " "
    #     + species_info["species"]
    #     + " "
    #     + species_info["strain"]
    # )
    source = pd.merge(
        isolation_src.loc[list(gp_locustag.columns), :],
        # pd.merge(
            
        #     species_info.loc[
        #         species_info.index.isin(list(gp_locustag.columns)), "genome_name"
        #     ],
        #     on=["genome_id"],
        # ),
        phylo_group,
        on=["genome_id"],
    ).sort_values("cluster")
    gp_locustag_ordered = gp_locustag.loc[:, list(source.index)]
    gp_binary_ordered = gp_binary.loc[list(source.index), :]

    n_genome = gp_binary.shape[0]

    logging.info(f"Number of genome: {n_genome}")
    summary["frequency"] = summary["No. isolates"] / n_genome

    test = pd.merge(
        gp_binary.loc[list(source.index), :], source.cluster, on=["genome_id"]
    )
    test.columns.names = ["Gene"]
    test_result = pd.merge(
        test.groupby(["cluster"]).mean().T, summary["pangenome_class_2"], on=["Gene"]
    )
    test_result_final = test_result.loc[
        test_result["pangenome_class_2"] == "Accessory",
        test_result.columns != "pangenome_class_2",
    ]

    # filter out the intra-species analysis genes
    final = test_result_final.apply(filter_cluster, axis=1)
    target_genes = final[final != False]
    target_gene_ordered = target_genes.sort_values().index

    # Target genes heatmap
    result_dict_target = {
        "rows": rows(source),
        "cols": cols(
            summary.loc[target_gene_ordered, :],
            annotation.loc[target_gene_ordered, :],
            gp_locustag_ordered.loc[target_gene_ordered, :],
        ),
        "matrix": matrix(gp_binary_ordered.loc[:, target_gene_ordered]),
    }

    with gzip.open(f"{heatmap_base_path}target.json.gz", "wt") as handle:
        json.dump(result_dict_target, handle)

    # gene_count = 0

    for gene_class in ['core', 'accessory', '1_15', 'above_1', 'only_1']:
        if gene_class == 'core':
            freq_filtered_gene = list(summary.loc[summary['pangenome_class_2'] == 'Core', :].index)
        elif gene_class == 'accessory':
            freq_filtered_gene = list(summary.loc[summary['pangenome_class_2'] == 'Accessory', :].index)
        elif gene_class == '1_15':
            freq_filtered_gene = list(summary.loc[(summary['pangenome_class_2'] == 'Rare') & (summary['frequency'] >= 0.01), :].index)
        elif gene_class == 'above_1':
            freq_filtered_gene = list(summary.loc[(summary['pangenome_class_2'] == 'Rare') & (summary['frequency'] < 0.01) & (summary['No. isolates'] > 1), :].index)
        else:
            freq_filtered_gene = list(summary.loc[summary['No. isolates'] == 1, :].index)

        with open(f"{sourceinfo_base_path}{gene_class}.json", 'w') as f:
            json.dump(source_info(source), f)

        result_dict_target = {"rows":rows(source),
                        "cols":cols(summary.loc[freq_filtered_gene,:],
                                    annotation.loc[freq_filtered_gene,:],
                                    gp_locustag_ordered.loc[freq_filtered_gene,:]),
                        "matrix":matrix(gp_binary_ordered.loc[:,freq_filtered_gene])}
        with gzip.open(f"{heatmap_base_path}{gene_class}.json.gz", "wt") as handle:
            json.dump(result_dict_target, handle)

    #     with open('../web_data/species/' + species + '/heatmap' + '_' + gene_class + '.json', 'w') as file:
    #         json.dump(result_dict_target, file)

    #     gzip_file('../web_data/species/' + species + '/heatmap' + '_' + gene_class + '.json', '../web_data/species/' + species + '/heatmap' + '_' + gene_class + '.json.gz')

    #     gene_count += len(freq_filtered_gene)


def run(args):
    generate_heatmap(
        args.gp_binary,
        args.gp_locustag,
        args.summary,
        args.eggnog_summary,
        args.mash_list,
        args.isosource,
        # args.species_info,
        args.heatmap_base,
        args.sourceinfo_base,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
