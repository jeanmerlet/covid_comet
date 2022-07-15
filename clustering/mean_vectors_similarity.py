import os
import pandas as pd
import numpy as np

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl/mean_vectors'
vecs_path = os.path.join(data_dir, 'mean_vectors-1.2')
out_path = os.path.join(data_dir, 'mean_vectors_edges.tsv')

vecs = pd.read_csv(vecs_path, sep='\t', header=None, index_col=None)

num_clusters = vecs.shape[0]
for i in range(num_clusters):
    for j in range(num_clusters):
        if j >= i: continue
        vec_a = vecs.iloc[i, :].values
        vec_b = vecs.iloc[j, :].values
        pearson = np.corrcoef(vec_a, vec_b)[0, 1]
        with open(out_path, 'a') as out_file:
            out_file.write(str(i+1) + '\t' + str(j+1) + '\t' + str(round(pearson, 2)) + '\n')

