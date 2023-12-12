state_tped_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/cali/ns-as-0s_mutation_count_filtered_100/combined_cali.tped'

out_tped_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/cali/ns-as-0s_mutation_count_filtered_100/combined_cali_counties-only.tped'

valid_seqs_path = '/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/metadata/metadata_2022_06_02/cali-counties/cali-counties_epi-ids.txt'

valid_seqs = {}
with open(valid_seqs_path, 'r') as in_file:
    for line in in_file:
        seq_id = line.strip()[8:].zfill(9)
        valid_seqs[seq_id] = True

with open(out_tped_path, 'a') as out_file:
    with open(state_tped_path, 'r') as in_file:
        for line in in_file:
            split_line = line.strip().split('\t')
            seq_id = split_line[1]
            if valid_seqs.get(seq_id):
                out_file.write(line)
