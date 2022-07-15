import numpy as np
import pandas as pd
from mpi4py import MPI
import os, re, time
from itertools import compress

root_in_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/transposed'
root_out_dir = os.path.join(root_in_dir, 'ns-as-0s_mutation_count_filtered_100')
counts_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/mutation_counts/combined_mutation_counts.tsv'
try:
    os.mkdir(root_out_dir)
except OSError:
    pass

nt_conversions = {'a': 'a', 'c': 'c', 'g': 'g', 't': 't', '-': '-', 'n': 'n',
                  'u': 'n', 'r': 'n', 'y': 'n', 'k': 'n', 'm': 'n', 's': 'n',
                  'w': 'n', 'b': 'n', 'd': 'n', 'h': 'n', 'v': 'n'}
nt_conversions = pd.Series(nt_conversions)

encoding_scheme = {'n': 'A\tA\tA\tA\tA',
                   'a': 'A\tA\tA\tA\tT',
                   't': 'A\tA\tA\tT\tA',
                   'c': 'A\tA\tT\tA\tA',
                   'g': 'A\tT\tA\tA\tA',
                   '-': 'T\tA\tA\tA\tA'}
encoding_scheme = pd.Series(encoding_scheme)

mutation_counts = pd.read_csv(counts_path, sep='\t', index_col=0, header=None)
pass_idx = mutation_counts.loc['valid', :].values == 'True'
pass_idx = list(compress(list(range(len(pass_idx))), pass_idx))

tsv_names = []
for r, d, f in os.walk(root_in_dir):
    for tsv in f:
        if re.search('^transposed_preprocessed.*tsv$', tsv):
            tsv_names.append(os.path.join(r, tsv))

tsv_names.sort()
num_sets = len(tsv_names)

def check_mutations(seq, mut_counts, align_seq):
    mutation_idx = seq != align_seq
    mutation_idx = np.arange(len(mutation_idx))[mutation_idx]
    for idx in mutation_idx:
        nt = seq[idx]
        if nt != 'n' and mut_counts.loc[nt, idx] < 100:
            seq[idx] = align_seq[idx]
    return seq

def convert_bases_to_tped_bin(seq_bases, encoding_scheme):
    seq_bases = encoding_scheme[seq_bases].values
    seq_bases = '\t'.join(seq_bases)
    return seq_bases

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

start = time.time()
for i, seq_set in enumerate(tsv_names):
    # distribute across ranks
    if i % size != rank: continue
    if i == 0: start_time = time.time()
    _, seq_set_name = os.path.split(seq_set)
    out_name = os.path.join(root_out_dir, 'mutation-filtered_' + seq_set_name[:-4] + '.tped')
    # don't redo a file if it's already done
    if os.path.isfile(out_name): continue
    print(f'filtering {seq_set_name} ({i+1}/{num_sets})')

    with open(seq_set) as in_file:
        for line_num, line in enumerate(in_file):
            if line_num == 0: continue
            if line_num - 1 not in pass_idx: continue
            nts = line.strip().lower().split('\t')
            nt_pos, nts = nts[0], np.array(nts[1:])
            nts = nt_conversions[nts].values
            if i == 0:
                seq_header = '0\t' + str(nt_pos) + '\t0\t0\t'
            else:
                seq_header = ''
            nts = convert_bases_to_tped_bin(nts, encoding_scheme)
            with open(out_name, 'a') as out_file:
                out_file.write(seq_header + nts + '\n')

    if i == 0: print(f'\ntotal time: {time.time() - start_time}')
