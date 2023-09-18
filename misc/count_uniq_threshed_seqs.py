# run this on the catted HipMCL input file (which is the threshed CoMet output)

edge_file = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/thresh_run_HH-0.19991/out_combined.tsv'

out_just_ids_file = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/thresh_run_HH-0.19991/all_ids.txt'

with open(out_just_ids_file, 'wt') as out_file:
    with open(edge_file, 'rt') as in_file:
        for i, line in enumerate(in_file):
            seq_id1, seq_id2, _ = line.strip().split('\t')
            seq_id1 = seq_id1[:-2]
            seq_id2 = seq_id2[:-2]
            out_file.write(seq_id1 + '\n')
            out_file.write(seq_id2 + '\n')

print(f'total number of lines (edgtes) in input file: {i}')
