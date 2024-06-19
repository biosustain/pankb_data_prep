import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import json
from pathlib import Path

def remove_slash(s):
    if '/' in str(s):
        return(s.replace('/','_'))
    else:
        return(s)

def generate_locustag_data(gp_locustag_path, all_locustag_path, gene_locustag_dir):
    df_gene_presence_locustag = pd.read_csv(gp_locustag_path, index_col='Gene',low_memory=False)
    df_gene_presence_locustag.index = [remove_slash(i) for i in list(df_gene_presence_locustag.index)]
    all_locustag_df = pd.read_csv(all_locustag_path, index_col=0, low_memory=False)
    
    for gene_id in df_gene_presence_locustag.index.tolist():    
        gene_locustag = []
        for genome_id, locus_tag in (df_gene_presence_locustag.loc[gene_id, :].dropna()).items():
            if "\t" in locus_tag:
                locus_tag_list = locus_tag.split("\t")
                for locus_tag_id in locus_tag_list:
                    gene_locustag.append(f"{genome_id}@{locus_tag_id}")
            else:
                gene_locustag.append(f"{genome_id}@{locus_tag}")

        # Write the JSON object to a file
        gene_locustag_dir = Path(gene_locustag_dir)
        with open(gene_locustag_dir / f"{gene_id}.json", 'w') as f:
            json.dump(all_locustag_df.loc[gene_locustag, :].to_dict(orient="records"), f, separators=(',', ':'), ensure_ascii=False)