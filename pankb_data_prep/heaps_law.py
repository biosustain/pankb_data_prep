import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error
import math
from scipy.stats import linregress
import json


def heaps_law(gp_binary_path, gene_freq_path):
    gp_binary = pd.read_csv(gp_binary_path, index_col=0)

    # -------------------------------------------------------------------------------------------------------------- #

    # Run the simulation

    # Sort the matrix by the sum of strains presence
    df = pd.DataFrame()
    df["data"] = gp_binary.sum(axis=1).sort_values(ascending=True)

    # Group the data by value and count the number of occurrences
    df = df.groupby('data').size().reset_index(name='frequency')

    # Sort the data in ascending order
    df = df.sort_values('data')

    # Calculate the cumulative frequency
    df['cumulative_frequency'] = df['frequency'].cumsum()

    # Calculate the core, accessory, and pangenome sizes
    x15 = math.floor(df['data'].quantile(0.15))
    x99 = math.floor(df['data'].quantile(0.99).round()) - 1
    core_size = df.loc[df['data'] >= x99, 'frequency'].sum()
    accessory_size = df[(df['data'] < x99) & (df['data'] >= x15)]['frequency'].sum() 
    rare_size = df.loc[df['data'] < x15, 'frequency'].sum()

    # Add Heap's law plot
    # Example gene presence matrix
    gene_presence_matrix = gp_binary.T
    # Set the number of genomes
    num_genomes = gene_presence_matrix.shape[0]

    # Create a list to store the Heap's law curves for each simulation
    pan_curves = []
    core_curves = []
    acc_curves = []

    # Set the number of simulations
    num_sims = 30

    # Run the simulations
    for i in range(num_sims):
        # Shuffle the rows of the gene presence matrix
        gene_presence_matrix_shuffled = gene_presence_matrix.sample(frac=1, replace=False)

        # Calculate Heap's law curve for pangenome
        x_pan = np.arange(num_genomes) + 1
        y_pan = np.zeros(num_genomes)
        for j in range(num_genomes):
            y_pan[j] = np.sum(np.any(gene_presence_matrix_shuffled.iloc[:j+1], axis=0))

        # Calculate Heap's law curve for core genome
        x_core = np.arange(num_genomes) + 1
        y_core = np.zeros(num_genomes)
        for j in range(num_genomes):
            y_core[j] = np.sum(np.all(gene_presence_matrix_shuffled.iloc[:j+1], axis=0))


        # Calculate Heap's law curve for accessory genome
        x_acc = np.arange(num_genomes) + 1
        y_acc = np.zeros(num_genomes)
        for j in range(num_genomes):
            y_acc[j] = np.sum((np.sum(gene_presence_matrix_shuffled.iloc[:j+1], axis=0) > np.ceil((j+1) * 0.15)) &
                            (np.sum(gene_presence_matrix_shuffled.iloc[:j+1], axis=0) <= np.floor((j+1) * 0.99))) + y_core[j]   

        # Append the curves to the lists
        pan_curves.append(y_pan)
        core_curves.append(y_core)
        acc_curves.append(y_acc)

    # Calculate the average Heap's law curves
    avg_pan = np.mean(np.array(pan_curves), axis=0)
    avg_core = np.mean(np.array(core_curves), axis=0)
    avg_acc = np.mean(np.array(acc_curves), axis=0)

    # -------------------------------------------------------------------------------------------------------------- #

    # Assuming you have the x and y values for the Heap's law curve
    x = np.arange(num_genomes) + 1  # Array of genome numbers
    y = avg_pan  # Array of unique gene counts

    # Perform linear regression
    slope, _, _, _, _ = linregress(np.log(x), np.log(y))

    # The slope represents the exponent in Heap's law equation (β)
    # logging.info("Slope (β):", slope)

    # -------------------------------------------------------------------------------------------------------------- #

    # Save data
    df.index = df.data
    data_dict = {'x_core':x_core.tolist(), 'avg_core':avg_core.tolist(),'x_acc':x_acc.tolist(),
                'avg_acc':avg_acc.tolist(),'x_pan':x_pan.tolist(),'avg_pan':avg_pan.tolist(),
                'frequency':(gp_binary.sum(axis=1).sort_values(ascending=False)).tolist(),'x15':int(x15),'x99':int(x99),
                'core_size':int(core_size),'accessory_size':int(accessory_size),'rare_size':int(rare_size),
                'rare_max_freq':int(df.loc[0:x15,'frequency'].max()),'accessory_max_freq':int(df.loc[x15:x99,'frequency'].max()), 
                'core_max_freq':int(df.loc[x99:,'frequency'].max()),'max_frequency_count':int(df['frequency'].max()),
                'data':(df['data']).tolist(),'cumulative_frequency':(df['cumulative_frequency']).tolist(),
                'max_cumulative':int(df['cumulative_frequency'].max()),'rare_max_cumulative':int(df.loc[:x15,'cumulative_frequency'].max()),
                'accessory_max_cumulative':int(df.loc[x15:x99,'cumulative_frequency'].max()),
                'core_max_cumulative':int(df.loc[x99:,'cumulative_frequency'].max()),
                'alpha': 1-slope}

    with open(gene_freq_path, 'w') as f:
        json.dump(data_dict, f)
        