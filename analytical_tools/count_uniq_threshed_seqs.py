# run this on the catted HipMCL input file (which is the threshed CoMet output)

edge_file = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/thresh_run_HH-0.19991/out_combined.tsv'

uniq_ids = set()

with open(edge_file) as in_file:
    for i, line in enumerate(in_file):
        seq_id1, seq_id2, _ = line.strip().split('\t')
        seq_id1 = seq_id1[:-2]
        seq_id2 = seq_id2[:-2]
        uniq_ids.add(seq_id1)
        uniq_ids.add(seq_id2)
        if i % 10000000000 == 0: print(i)

print(len(uniq_ids))
