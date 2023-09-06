import pandas as pd
import numpy as np
import re
import os

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/mutation_counts'

mutation_cutoff = 10000

tsv_names = []
for r, d, f in os.walk(data_dir):
    for tsv in f:
        if 'combined' not in tsv:
            tsv_names.append(os.path.join(r, tsv))

total_counts = pd.read_csv(tsv_names[0], index_col=0, header=None, sep='\t')
tsv_names = tsv_names[1:]

for tsv_name in tsv_names:
    counts = pd.read_csv(tsv_name, index_col=0, header=None, sep='\t')
    total_counts += counts

total_counts.columns = np.arange(total_counts.shape[1])
total_counts.loc['valid', :] = (total_counts >= mutation_cutoff).any(axis=0)

out_path = os.path.join(data_dir, f'combined_mutation_counts_{mutation_cutoff}.tsv')
total_counts.to_csv(out_path, sep='\t', header=None)
