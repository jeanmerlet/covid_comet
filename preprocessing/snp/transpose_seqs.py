import os, time
import numpy as np
import pandas as pd
from mpi4py import MPI

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665'
out_dir = os.path.join(data_dir, 'transposed')

try:
    os.mkdir(out_dir)
except OSError:
    pass

seq_paths = []
for seq_name in os.listdir(data_dir):
    if '.tsv' in seq_name:
        seq_paths.append(os.path.join(data_dir, seq_name))

def write_transpose(in_path, out_path, remove_ref):
    data = pd.read_csv(in_path, sep='\t', header=None, index_col=None)
    if remove_ref:
        data = data.drop(axis=0, index=0)
    data = data.T
    data.columns = data.iloc[0, :].values
    data = data.drop(axis=0, index=0)
    data.index = list(range(0, data.shape[0]))
    data.to_csv(out_path, sep='\t', index=True, header=True)
    del(data)
    
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

num_paths = len(seq_paths)
for i, seq_path in enumerate(seq_paths):
    # distribute across ranks
    if i % size != rank: continue
    _, seq_name = os.path.split(seq_path)
    print(f'transposing {seq_name} ({i+1}/{num_paths})\n')
    out_path = os.path.join(out_dir, 'transposed_' + seq_name)
    remove_ref = False
    if '001' not in seq_name: remove_ref = True
    write_transpose(seq_path, out_path, remove_ref)
