import pandas as pd
import numpy as np
from scipy.stats import linregress

def get_genome_list(args):
    with open(args.genomes, "r") as f:
        l = [genome.strip() for genome in f.readlines()]
    return l


COG_TABLE = pd.DataFrame.from_dict(
    {
        "A": "RNA processing and modification",
        "B": "Chromatin Structure and dynamics",
        "C": "Energy production and conversion",
        "D": "Cell cycle control and mitosis",
        "E": "Amino Acid metabolis and transport",
        "F": "Nucleotide metabolism and transport",
        "G": "Carbohydrate metabolism and transport",
        "H": "Coenzyme metabolis",
        "I": "Lipid metabolism",
        "J": "Tranlsation",
        "K": "Transcription",
        "L": "Replication and repair",
        "M": "Cell wall/membrane/envelop biogenesis",
        "N": "Cell motility",
        "O": "Post-translational modification, protein turnover, chaperone functions",
        "P": "Inorganic ion transport and metabolism",
        "Q": "Secondary Structure",
        "T": "Signal Transduction",
        "U": "Intracellular trafficing and secretion",
        "Y": "Nuclear structure",
        "Z": "Cytoskeleton",
        "R": "General Functional Prediction only",
        "S": "Function Unknown",
        "-": "Not found in COG",
    },
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