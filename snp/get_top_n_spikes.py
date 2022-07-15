import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os, re

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/transposed/binned'

count_paths = [os.path.join(data_dir, path) for path in os.listdir(data_dir)]
count_paths.sort()

for i, count_path in enumerate(count_paths):
    if i == 0:
        count_sums = pd.read_csv(count_path, sep='\t', index_col=0, header=0)
    else:
        counts = pd.read_csv(count_path, sep='\t', index_col=0, header=0)
        count_sums += counts

count_sums = count_sums.loc[:, ['-', 'n']]

def get_n_spikes(count_sums, num_spikes):
    count_sums = count_sums.sort_values(by='n', ascending=False)
    print(count_sums.iloc[:num_spikes, :])


num_spikes = 20
get_n_spikes(count_sums, num_spikes)
