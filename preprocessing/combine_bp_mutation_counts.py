import pandas as pd
import numpy as np
import re, os
import argparse


# arpgarse
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir', type=str, required=True)
parser.add_argument('-m', '--mutation_threshold', type=int, required=False, default=100)

args = parser.parse_args()
data_dir = args.data_dir
mutation_threshold = args.mutation_threshold


# mutation count paths
tsv_names = []
for r, d, f in os.walk(data_dir):
    for tsv in f:
        if 'combined' not in tsv:
            tsv_names.append(os.path.join(r, tsv))

total_counts = pd.read_csv(tsv_names[0], index_col=0, header=None, sep='\t')
tsv_names = tsv_names[1:]

# combine and write
for tsv_name in tsv_names:
    counts = pd.read_csv(tsv_name, index_col=0, header=None, sep='\t')
    total_counts += counts

total_counts.columns = np.arange(total_counts.shape[1])
total_counts.loc['valid', :] = (total_counts >= mutation_threshold).any(axis=0)

out_path = os.path.join(data_dir, f'combined_mutation_counts_{mutation_threshold}.tsv')
total_counts.to_csv(out_path, sep='\t', header=None)
