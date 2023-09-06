#!/bin/bash

in_txt_path="/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/2-way-0.965_mut-100/combined_2-way-0.965.txt"

out_tsv_path="/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/2-way-0.965_mut-100/combined_2-way-0.965.tsv"

cat $in_txt_path | tr ' ' '\t' | cut -f 5-7 > $out_tsv_path
