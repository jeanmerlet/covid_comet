import numpy as np

clusters_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/cali/ns-as-0s_mutation_count_filtered_100/comet_out/hip_mcl/clusters_inflation-1.2.tsv'

tped_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/cali/ns-as-0s_mutation_count_filtered_100/combined_cali.tped'

out_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/cali/ns-as-0s_mutation_count_filtered_100/comet_out/hip_mcl/out_seqs_mtx.tsv'

seq_ids = {}
with open(clusters_path, 'r') as in_file:
    for line in in_file:
        cluster = line.strip().split(' ')
        for seq_id in cluster:
            seq_ids[seq_id[:-2]] = True

print(len(seq_ids.keys()))
raise SystemExit()

with open(out_path, 'a') as out_file:
    with open(tped_path, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split('\t')
            seq_id = split_line[1]
            if seq_ids.get(seq_id) is not None:
                out_line = '\t'.join([seq_id, '\t'.join(split_line[4:])])
                out_file.write(out_line)
