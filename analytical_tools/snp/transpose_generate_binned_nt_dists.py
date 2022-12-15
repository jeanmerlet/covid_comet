import os, re, time
import numpy as np
import pandas as pd
from mpi4py import MPI

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/transposed'
out_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/transposed/binned'

snp_set_paths = []
for path in os.listdir(data_dir):
    if '.tsv' in path:
        snp_set_paths.append(os.path.join(data_dir, path))

snp_set_paths.sort()
num_sets = len(snp_set_paths)

codes = {'a': 0, 'c': 0, 'g': 0, 't': 0,
         'u': 0, 'r': 0, 'y': 0, 'k': 0,
         'm': 0, 's': 0, 'w': 0, 'b': 0,
         'd': 0, 'h': 0, 'v': 0, 'n': 0,
         '-': 0}

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, snp_set_path in enumerate(snp_set_paths):
    if i == 0: start = time.time()
    # distribute across ranks
    if i % size != rank: continue
    _, snp_set_name = os.path.split(snp_set_path)
    out_path = os.path.join(out_dir, 'nt-counts_' + snp_set_name)
    print(f'\ncounting nts in {snp_set_name} ({i+1}/{num_sets})')
    num_pos = 29324
    code_bins = [codes.copy() for i in range(num_pos)]

    with open(snp_set_path) as in_file:
        for i, line in enumerate(in_file):
            if i > 0:
                nts = line.strip().split('\t')[1:]
                nt_ids, counts = np.unique(nts, return_counts=True)
                for j, nt_id in enumerate(nt_ids):
                    code_bins[i-1][nt_id] += counts[j]

    out_df = pd.DataFrame(code_bins)
    out_df.to_csv(out_path, sep='\t')
    if i == 0: print(f'\n total time: {time.time() - start}')
