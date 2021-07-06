import os
import re
import numpy as np
import pandas as pd
from mpi4py import MPI
import time
import argparse

codes = {'a': 0, 'c': 0, 'g': 0, 't': 0,
         'u': 0, 'r': 0, 'y': 0, 'k': 0,
         'm': 0, 's': 0, 'w': 0, 'b': 0,
         'd': 0, 'h': 0, 'v': 0, 'n': 0,
         '-': 0}

gpfs_root = '/gpfs/alpine/syb105/proj-shared'

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', default=os.path.join(gpfs_root, 'Projects/GeoBio_CoMet/data/original'),
                    help='path to dataset')
parser.add_argument('-p', '--hosts', default=os.path.join(gpfs_root,
                    'Personal/jmerlet/projects/sars_cov_2_geo/data/metadata/non_human_hosts.tsv'),
                    help='path to non-human hosts metadata')
parser.add_argument('-g', '--d_cutoff', type=int, default=1000,
                    help='maximum number of indeterminate gaps in a sequence')
parser.add_argument('-n', '--n_cutoff', type=float, default=0.01,
                    help='maximum proportion of ns in a sequence')
parser.add_argument('-l', '--left_trim', type=int, default=342,
                    help="trim all bp locations below this")
parser.add_argument('-r', '--right_trim', type=int, default=29665,
                    help="trim all bp locations above this")
args = parser.parse_args()

out_dir_name = f'preprocessed_d-cutoff_{args.d_cutoff}_n-cutoff_{args.n_cutoff}_pos_{args.left_trim}-{args.right_trim}'
out_dir = os.path.join(args.data, out_dir_name)
log_dir = os.path.join(out_dir, 'logs')
mut_dir = os.path.join(out_dir, 'mutation_counts')
try:
    os.mkdir(out_dir)
except OSError:
    pass
try:
    os.mkdir(log_dir)
except OSError:
    pass
try:
    os.mkdir(mut_dir)
except OSError:
    pass

non_human = pd.read_csv(args.hosts, sep='\t', header=0, index_col=0).index.values

fasta_names = []
for r, d, f in os.walk(args.data):
    for fasta in f:
        if re.search('.*fasta$', fasta):
            fasta_names.append(os.path.join(r, fasta))

fasta_names.sort()
num_sets = len(fasta_names)

def pass_filters(seq_name, log_path, nt_counts, non_human,
                 wrong_code, d_cutoff, n_cutoff):
    d_prop = nt_counts['-']
    n_prop = nt_counts['n'] / seq_length
    fail_line = None
    if wrong_code:
        fail_line = f'"{seq_name}" contains an INVALID NT CODE: {wrong_code}'
    elif d_prop > d_cutoff:
        fail_line = f'"{seq_name}" FAILED to pass DELETION filter'
    elif n_prop > n_cutoff:
        fail_line = f'"{seq_name}" FAILED to pass N filter'
    elif re.search('(EPI_ISL_\d+)', seq_name) and re.search('(EPI_ISL_\d+)', seq_name).groups()[0] in non_human:
        fail_line = f'"{seq_name}" has NON-HUMAN HOST'

    if fail_line:
        with open(log_path, 'a') as out_log:
            out_log.write(fail_line + '\n')
        return False
    return True

def count_bases(nt_counts, seq, seq_length, wrong_code):
    bases = list(seq)
    seq_length += len(bases)
    uniq_bases, counts = np.unique(bases, return_counts=True)
    for i, nt_id in enumerate(uniq_bases):
        try:
            nt_counts[nt_id] += counts[i]
        except KeyError as key:
            wrong_code = key

    return nt_counts, seq_length, wrong_code

def preprocess(seq_name, sequence, out_path, left_trim_pos, right_trim_pos):
    sequence = list(''.join(sequence))
    sequence = sequence[left_trim_pos:right_trim_pos + 1]
    out_sequence = seq_name + '\t' + '\t'.join(sequence) + '\n'
    with open(out_path, 'a') as out_file:
        out_file.write(out_sequence)
    return(sequence)

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

for i, seq_set in enumerate(fasta_names):
    # distribute across ranks
    if i % size != rank: continue
    _, seq_set_name = os.path.split(seq_set)
    # time the largest file
    if re.search('04-12-2021', seq_set_name): total_start = time.time()
    out_path = os.path.join(out_dir, 'preprocessed_' + seq_set_name[:-6] + '.tsv')
    log_path = os.path.join(log_dir, seq_set_name[:-6] + '_preprocessing_log.txt')
    mutation_counts_out_path = os.path.join(mut_dir, seq_set_name[:-6] + '_mutation_counts.tsv')
    # don't redo a file if it's already done
    if os.path.isfile(out_path): continue
    print(f'filtering {seq_set_name} ({i+1}/{num_sets})\n')

    line_num = 0
    first_seq = True
    with open(seq_set) as in_file:
        for line in in_file:
            line_num += 1
            # sequence headers start with >
            if '>' in line:
                # check filters for an entire sequence at the start of the next sequence
                # preprocess and write the sequence if it passes all filters
                if line_num > 1:
                    if pass_filters(seq_name, log_path, nt_counts, non_human,
                                    wrong_code, args.d_cutoff, args.n_cutoff):
                        sequence = preprocess(seq_name, sequence, out_path, args.left_trim, args.right_trim)
                        if first_seq:
                            # first sequence is alignment seq
                            # and is stored for comparison
                            length = len(sequence)
                            align_seq = sequence.copy()
                            # also create array to count mutations
                            row_names = ['a', 'c', 'g', 't', '-']
                            data = np.zeros((len(row_names), length), dtype=np.int)
                            diffs = pd.DataFrame(data=data, index=row_names, columns=None)
                            first_seq = False
                        else:
                            diffs = count_seq_differences(diffs, row_names, align_seq, sequence, length)
                seq_name = line.strip()
                nt_counts = codes.copy()
                sequence = []
                seq_length = 0
                wrong_code = False
                continue
            else:
                # some sequences are uppercase
                line = line.strip().lower()
                sequence.append(line)
                nt_counts, seq_length, wrong_code = count_bases(nt_counts, line, seq_length, wrong_code)

    diffs.to_csv(mutation_counts_out_path, sep='\t', header=False)
    if re.search('04-12-2021', seq_set_name): print(f'total time: {time.time() - total_start}')
