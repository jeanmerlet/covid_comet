import os
import re
import numpy as np
import pandas as pd
from mpi4py import MPI
import time

root_data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original'
root_in_dir = os.path.join(root_data_dir, 'preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665')
root_out_dir = os.path.join(root_in_dir, 'ns-as-1s_mutation_count_filtered_100')
counts_path = os.path.join(root_in_dir, 'mutation_counts/combined_mutation_counts.tsv')
try:
    os.mkdir(root_out_dir)
except OSError:
    pass

nt_conversions = {'a': 'a', 'c': 'c', 'g': 'g', 't': 't', '-': '-', 'n': 'n',
                  'u': 'n', 'r': 'n', 'y': 'n', 'k': 'n', 'm': 'n', 's': 'n',
                  'w': 'n', 'b': 'n', 'd': 'n', 'h': 'n', 'v': 'n'}
nt_conversions = pd.Series(nt_conversions)

encoding_scheme = {'n': 'T\tT\tT\tT\tT',
                   'a': 'A\tA\tA\tA\tT',
                   't': 'A\tA\tA\tT\tA',
                   'c': 'A\tA\tT\tA\tA',
                   'g': 'A\tT\tA\tA\tA',
                   '-': 'T\tA\tA\tA\tA'}
encoding_scheme = pd.Series(encoding_scheme)

mutation_counts = pd.read_csv(counts_path, sep='\t', index_col=0, header=None)
pass_idx = mutation_counts.loc['valid', :].values == 'True'
passed_mut_counts = mutation_counts.loc[:, pass_idx]
passed_mut_counts = passed_mut_counts.iloc[:-1, :].astype('float').astype(np.int32)
passed_mut_counts.columns = np.arange(len(passed_mut_counts.columns))

tsv_names = []
for r, d, f in os.walk(root_in_dir):
    for tsv in f:
        if re.search('^preprocessed.*tsv$', tsv):
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

def create_seq_header(seq_id):
    epi_id = re.search('\|.+_(\d+)\|.*$', seq_id).groups()[0]
    id_length = len(epi_id)
    if id_length < 9:
        adj_zeros = '0' * (9 - id_length)
        epi_id = adj_zeros + epi_id
    try:    
        assert(len(epi_id) == 9)
    except AssertionError as e:
        print(f'ERROR: {epi_id} of incorrect length ({e} - should be 9)')
        return None
    seq_header = '0\t' + epi_id + '\t0\t0'
    return seq_header

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
    _, seq_set_name = os.path.split(seq_set)
    out_name = os.path.join(root_out_dir, 'mutation-filtered_' + seq_set_name[:-4] + '.tped')
    # don't redo a file if it's already done
    if os.path.isfile(out_name): continue
    print(f'filtering {seq_set_name} ({i+1}/{num_sets})')
    # record time for largest file
    if '04-12-2021' in seq_set_name: start_time = time.time()

    first_seq = True
    with open(seq_set) as in_file:
        for line in in_file:
            seq = line.strip().lower().split('\t')
            seq_id, seq = seq[0], np.array(seq[1:])
            seq = seq[pass_idx]
            seq_header = None
            if first_seq:
                first_seq = False
                # check it's the alignment sequence
                assert '_045512.2' in seq_id
                align_seq = seq.copy()
                # wuhan sequence is duplicated in each file
                # and doesn't have the same header pattern
                if i == 0:
                    seq_header = '0\t000000000\t0\t0'
                else:
                    continue
            seq = nt_conversions[seq].values
            seq = check_mutations(seq, passed_mut_counts, align_seq)
            if not seq_header:
                seq_header = create_seq_header(seq_id)
            if seq_header:
                seq = convert_bases_to_tped_bin(seq, encoding_scheme)
                seq = seq_header + '\t' + seq + '\n'
                with open(out_name, 'a') as out_file:
                    out_file.write(seq)

    if '04-12-2021' in seq_set_name: print(f'\ntotal time: {time.time() - start_time}')
