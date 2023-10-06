import os
import re
import numpy as np
import pandas as pd
from mpi4py import MPI
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--root_in_dir', type=str, required=True)
parser.add_argument('-n', '--name', type=str, required=True)
parser.add_argument('-w', '--seq_id_whitelist_path', type=str, required=True)
parser.add_argument('-p', '--file_prefix', type=str, required=True)

args = parser.parse_args()
root_in_dir = args.root_in_dir
name = args.name
seq_id_whitelist_path = args.seq_id_whitelist_path
file_prefix = args.file_prefix


root_out_dir = os.path.join(root_in_dir, name)
mut_dir = os.path.join(root_out_dir, 'mutation_counts')
try:
    os.mkdir(root_out_dir)
except OSError:
    pass
try:
    os.mkdir(mut_dir)
except OSError:
    pass


seq_id_whitelist = np.squeeze(pd.read_csv(seq_id_whitelist_path).values)

tsv_names = []
for r, d, f in os.walk(root_in_dir):
    for tsv in f:
        if re.search('^preprocessed.*tsv$', tsv):
            tsv_names.append(os.path.join(r, tsv))
tsv_names.sort()
num_sets = len(tsv_names)


#remake mutation counts
def count_seq_differences(diffs, nt_codes_to_check, seq1, seq2, length):
    for i in range(length):
        nt_value = seq2[i]
        if nt_value in nt_codes_to_check:
            if seq1[i] != nt_value:
                diffs.loc[nt_value, i] += 1
    return diffs


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

start = time.time()
for i, seq_set in enumerate(tsv_names):
    # distribute across ranks
    if i % size != rank: continue
    _, seq_set_name = os.path.split(seq_set)
    out_name = os.path.join(root_out_dir, file_prefix + '_' + seq_set_name)
    mutation_counts_out_path = os.path.join(mut_dir, 'mutation_counts_' + seq_set_name)
    print(f'subsetting {seq_set_name} ({i+1}/{num_sets})', flush=True)
    # record time for first file
    if i == 0: start_time = time.time()
    first_seq = True
    with open(seq_set) as in_file:
        for line in in_file:
            seq = line.strip().split('\t')
            seq_id, seq = seq[0], seq[1:]
            # wuhan refseq is first seq and doesn't have an epi-id
            # write it every time then subset the other seqs
            if first_seq:
                first_seq = False
                with open(out_name, 'a') as out_file:
                    out_file.write(line)
                # also create array to count mutations
                row_names = ['a', 'c', 'g', 't', '-']
                length = len(seq)
                data = np.zeros((len(row_names), length), dtype=int)
                diffs = pd.DataFrame(data=data, index=row_names, columns=None)
                align_seq = seq.copy()
            # otherwise, write subsampled seqs
            else:
                seq_id = re.search('(EPI_ISL_\d+)', seq_id).groups(0)[0]
                # check if the epi-id is in the pre-made epi-id whitelist
                if seq_id not in seq_id_whitelist:
                    continue
                else:
                    with open(out_name, 'a') as out_file:
                        out_file.write(line)
                # for mutation counts
                diffs = count_seq_differences(diffs, row_names, align_seq, seq, length)
    diffs.to_csv(mutation_counts_out_path, sep='\t', header=False)
    if i == 0: print(f'\ntotal time: {time.time() - start_time}', flush=True)
