import pandas as pd
import numpy as np
from scipy.stats import linregress

def get_genome_list(args):
    with open(args.genomes, "r") as f:
        l = [genome.strip() for genome in f.readlines()]
    return l

COG_DICT = {
        "A": "RNA processing and modification",
        "B": "Chromatin structure and dynamics",
        "C": "Energy production and conversion",
        "D": "Cell cycle control and mitosis",
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
        "T": "Signal transduction",
        "U": "Intracellular trafficing and secretion",
        "V": "Defense mechanisms",
        "W": "Extracellular structures",
        "X": "Mobilome: prophages, transposons",
        "Y": "Nuclear structure",
        "Z": "Cytoskeleton",
        "R": "General function prediction only",
        "S": "Function unknown",
        "-": "Not found in COG",
    }
COG_TABLE = pd.DataFrame.from_dict(
    COG_DICT,
    orient="index",
    columns=["Function details"],
)
# COG_TABLE.index.name = "Categories"
COG_TABLE["Categories"] = COG_TABLE.index

def calculate_lambda(gene_presence_matrix, num_samples=30, num_repetitions=50):
    np.random.seed(1)
    num_genomes, num_genes = gene_presence_matrix.shape
    lambda_values = []

    for _ in range(num_repetitions):
        # Randomly sample 30 genomes
        sampled_matrix = gene_presence_matrix.sample(n=min(num_samples, num_genomes), replace=False)
        
        # Calculate pangenome size for increasing numbers of genomes
        x = np.arange(1, len(sampled_matrix) + 1)
        y = np.zeros(len(sampled_matrix))
        
        for i in range(len(sampled_matrix)):
            y[i] = np.sum(np.any(sampled_matrix.iloc[:i+1], axis=0))
        
        # Perform linear regression on log-transformed data
        slope, _, _, _, _ = linregress(np.log(x), np.log(y))
        
        # Calculate lambda 
        lambda_value = slope
        lambda_values.append(lambda_value)
    
    # Calculate mean and standard deviation of lambda values
    mean_lambda = np.mean(lambda_values)
    std_lambda = np.std(lambda_values)
    
    return mean_lambda, std_lambda