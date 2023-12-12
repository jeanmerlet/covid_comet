#!/bin/bash

in_txt_path="/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/cali/ns-as-0s_mutation_count_filtered_100/comet_out/combined_out.txt"

out_tsv_path="/lustre/orion/syb111/proj-shared/Projects/sars-cov-2_geo/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/cali/ns-as-0s_mutation_count_filtered_100/comet_out/combined_out.tsv"

cat $in_txt_path | tr ' ' '\t' | cut -f 5-7 > $out_tsv_path
