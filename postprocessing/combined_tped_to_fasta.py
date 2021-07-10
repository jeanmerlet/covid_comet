import pandas as pd
import numpy as np
import re

tped_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/mutation_count_filtered_100/small_combined/seq_from_each_40.tped'
seq_id_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/data/small_run_exploration/seq40_exported_hh_0.799_cyto_clusters.tsv'
out_prefix = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/mutation_count_filtered_100/small_combined/seq_from_each_40_mcl_cluster_fastas/seq_from_each_40.fasta'

decoding_scheme = {'A\tA\tA\tA\tA': 'n',
                   'A\tA\tA\tA\tT': 'a',
                   'A\tA\tA\tT\tA': 't',
                   'A\tA\tT\tA\tA': 'c',
                   'A\tT\tA\tA\tA': 'g',
                   'T\tA\tA\tA\tA': '-'}
decoding_scheme = pd.Series(decoding_scheme)

seq_ids = pd.read_csv(seq_id_path, sep='\t', header=0, index_col=None, dtype=str)

with open(tped_path) as in_file:
    for line in in_file:
        line = line.strip().split('\t')
        header, seq = line[:4], line[4:]
        seq_id = header[1]
        if seq_id in seq_ids['name'].values:
            new_seq = []
            seq = ['\t'.join(seq[i*5:i*5+5]) for i in range(len(seq) // 5)]
            seq = decoding_scheme[seq]
            cluster = seq_ids.loc[:, 'mcl_cluster'].values[seq_ids['name'] == seq_id][0]
            out_suffix = f'_cluster_{cluster}_hh_7.9'
            out_path = re.sub('\.fasta', out_suffix + '.fasta', out_prefix)
            with open(out_path, 'a') as out_file:
                out_file.write('>' + ''.join(header) + '\n')
                out_file.write(''.join(seq) + '\n')
