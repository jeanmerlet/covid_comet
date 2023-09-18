import os, re

out_text_files_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/thresh_run_HH-0.19991/txt'
out_path = os.path.join(out_text_files_dir, 'missing_files.txt')

all_files = os.listdir(out_text_files_dir)

with open(out_path, 'wt') as out_file:
    i = 0
    for f in all_files:
        if f[:3] == 'out':
            if 'out_' + str(i).zfill(4) + '.txt' not in all_files:
                out_file.write('out_' + str(i).zfill(4) + '.txt' + '\n')
            i += 1

      
