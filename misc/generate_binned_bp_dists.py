import os
import re
import numpy as np
import pandas as pd
from mpi4py import MPI
import time

root_data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/count_filtered_0.01N_1000D'
root_out_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/data/filtered_0.01N_1000D_counts'

try:
    os.mkdir(root_data_dir)
except OSError:
    pass

try:
    os.mkdir(root_out_dir)
except OSError:
    pass

codes = {'a': 0, 'c': 0, 'g': 0, 't': 0,
         'u': 0, 'r': 0, 'y': 0, 'k': 0,
         'm': 0, 's': 0, 'w': 0, 'b': 0,
         'd': 0, 'h': 0, 'v': 0, 'n': 0,
         '-': 0}

fasta_names = []
for r, d, f in os.walk(root_data_dir):
    for fasta in f:
        if re.search('.*fasta$', fasta):
            fasta_names.append(os.path.join(r, fasta))

num_sets = len(fasta_names)

# set bin size TODO: implement some kind of args system
bin_size = 10
n_bins = 30000 // bin_size
root_out_dir = os.path.join(root_out_dir, 'bin_size_' + str(bin_size))

try:
    os.mkdir(root_out_dir)
except OSError:
    pass

def bin_nt_counts(seq_name, seq, code_bins, bin_size):
    seq = ''.join(seq)
    num_bins = len(seq) // bin_size
    binned_seq = [seq[(i*bin_size):(i*bin_size + bin_size)] for i in range(num_bins)]
    for i, bases in enumerate(binned_seq):
        uniq_bases, counts = np.unique(list(bases), return_counts=True)
        for j, nt_id in enumerate(uniq_bases):
            code_bins[i][nt_id] += counts[j]

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, seq_set in enumerate(fasta_names):
    if i == 0: start = time.time()
    # distribute across ranks
    if i % size != rank: continue
    _, seq_set_name = os.path.split(seq_set)
    out_name = os.path.join(root_out_dir, seq_set_name[:-6] + '_nt-counts_' + str(bin_size) + '.tsv')
    print(f'\ncounting bases in {seq_set_name} ({i+1}/{num_sets})')
    code_bins = [codes.copy() for i in range(n_bins)]
    line_num = 0

    with open(seq_set) as in_file:
        for line in in_file:
            line_num += 1
            # sequences all start with the sequence id lines, which start with >
            if '>' in line:
                # bin and write the bins to file for the previous sequence
                if line_num > 1:
                    bin_nt_counts(seq_name, sequence, code_bins, bin_size)
                seq_name = line
                sequence = []
                continue
            else:
                # lowercase because some sequences are uppercase
                bases = line.strip().lower()
                sequence.append(bases)

    out_df = pd.DataFrame(code_bins)
    out_df.to_csv(out_name, sep='\t')
    if i == 0: print(f'\n total time: {time.time() - start}')
