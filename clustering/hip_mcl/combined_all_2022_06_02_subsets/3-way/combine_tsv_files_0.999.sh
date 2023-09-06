#!/bin/bash

in_tsv_dir="/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100/dedup"
out_tsv_path="/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100/dedup/combined_dedup_3-way.tsv"

cat $in_tsv_dir/*tsv > $out_tsv_path
