import os, re

in_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/combined.tsv'

out_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/combined_vite_format.tsv'

with open(in_path, 'r') as in_file:
    with open(out_path, 'w') as out_file:
        for i, line in enumerate(in_file):
            if i % 10000000 == 0: print(i, flush=True)
            seq_id1, seq_id2, value = line.strip().split('\t')
            seq_id1 = seq_id1[:-2]
            seq_id2 = seq_id2[:-2]
            if seq_id1 == '000000001' or seq_id2 == '000000001':
                print('FAIL', flush=True)
            if seq_id1 == '000000000':
                seq_id1 = '000000001'
            if seq_id2 == '000000000':
                seq_id2 = '000000001'
            out_file.write('\t'.join([seq_id1, seq_id2, value]) + '\n')
