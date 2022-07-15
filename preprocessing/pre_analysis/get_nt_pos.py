import numpy as np
import pandas as pd
import os

mut_counts_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/mutation_counts/combined_mutation_counts.tsv'
out_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2021_09_30/nt_pos.tsv'

counts = pd.read_csv(mut_counts_path, sep='\t', index_col=0, header=None)

nt_pos = counts.columns[counts.loc['valid', :].values == 'True'].values
nt_pos = list(nt_pos)
nt_pos = [str(x + 341) for x in nt_pos]
nt_pos = '\t'.join(nt_pos) + '\n'

with open(out_path, 'wt') as out_file:
    out_file.write(nt_pos)

