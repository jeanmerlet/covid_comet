import pandas as pd
import numpy as np
import re

tped_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/mutation_count_filtered_100/small_combined/seq_from_each_10.tped'

out_path = re.sub('tped', 'fasta', tped_path)

decoding_scheme = {'A\tA\tA\tA\tA': 'n',
                   'A\tA\tA\tA\tT': 'a',
                   'A\tA\tA\tT\tA': 't',
                   'A\tA\tT\tA\tA': 'c',
                   'A\tT\tA\tA\tA': 'g',
                   'T\tA\tA\tA\tA': '-'}
decoding_scheme = pd.Series(decoding_scheme)

with open(tped_path) as in_file:
    for line in in_file:
        line = line.strip().split('\t')
        header, seq = line[:4], line[4:]
        new_seq = []
        seq = ['\t'.join(seq[i*5:i*5+5]) for i in range(len(seq) // 5)]
        seq = decoding_scheme[seq]
        #line = '\t'.join(header) + '\t' + '\t'.join(seq) + '\n'
        with open(out_path, 'a') as out_file:
            out_file.write('>' + ''.join(header) + '\n')
            out_file.write(''.join(seq) + '\n')
