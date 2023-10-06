import os
import re
import pandas as pd
import numpy as np

tped_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100'
sequence_nums_path =  os.path.join(tped_dir, 'sequence_nums.tsv')
out_dir = os.path.join(tped_dir, 'small_combined')
try:
    os.mkdir(out_dir)
except OSError:
    pass
num_seqs_per_file = 40
out_path = os.path.join(out_dir, 'seq_from_each_' + str(num_seqs_per_file) + '.tped')

tped_paths = []
for r, d, f in os.walk(tped_dir):
    for tped in f:
        if re.search('^mutation.*tped', tped):
            tped_paths.append(os.path.join(r, tped))

tped_paths.sort()

sequence_nums = pd.read_csv(sequence_nums_path, header=None, index_col=0, sep='\t')

for i, tped_path in enumerate(tped_paths):
    print(i)
    head, tail = os.path.split(tped_path)
    num_seqs = np.arange(2, sequence_nums.loc[tail].values[0] + 1)
    seq_lines = np.random.choice(num_seqs, num_seqs_per_file, replace=False)
    with open(tped_path, 'r') as in_file:
        line_num = 0
        for line in in_file:
            line_num += 1
            # only get the alignment sequence from the first tped
            if i == 0 and line_num == 1:
                with open(out_path, 'a') as out_file:
                    out_file.write(line)
            elif line_num in seq_lines:
                with open(out_path, 'a') as out_file:
                    out_file.write(line)
