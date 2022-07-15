import os, sys

#input_path = sys.argv[1]
#threshold = float(sys.argv[2])

input_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/out_868.txt'
output_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/dedup'
threshold = 0
dup_counter = 0

_, input_filename = os.path.split(input_path)
with open(input_path + '.parsed', 'w') as output_file:
    with open(input_path) as input_file:
        edge_dict   = {}
        for line in input_file:
            _, _, _, _, seq_id1, seq_id2, weight = line.split(' ')
            #weight = float(line[6][:-1])
            if float(weight) >= threshold:
                if int(seq_id1[:-2]) > int(seq_id2[:-2]):
                    seq_id1, seq_id2 = seq_id2, seq_id1
                if (seq_id1, seq_id2) not in edge_dict:
                    edge_dict[(seq_id1, seq_id2)] = weight
                else:
                    dup_counter += 1
                    edge_dict[(seq_id1, seq_id2)] = max(weight, edge_dict[(seq_id1, seq_id2)])

    # write to file
    for edge in edge_dict:
        parsed_line = ' '.join((edge[0], edge[1], edge_dict[edge] + '\n'))
        output_file.write(parsed_line)

print(f'# of dups: {dup_counter}')
