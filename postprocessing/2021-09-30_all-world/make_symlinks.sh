in_dir=/gpfs/alpine/syb105/world-shared/joubert/tarsnew/combined_all-2022-06-03
out_dir=/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/frontier_run_out

for file in $(ls $in_dir/*bin)
do
    ln -s $file $out_dir
done
