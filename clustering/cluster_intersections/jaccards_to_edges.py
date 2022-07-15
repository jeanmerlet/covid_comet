import os, re, time
import pandas as pd
import numpy as np

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl/jaccard'
jaccard_path = os.path.join(data_dir, 'intersections_top_100_clusters.tsv')
out_path = os.path.join(data_dir, 'edges_top_100_clusters.tsv')

jaccard = pd.read_csv(jaccard_path, sep='\t', header=0, index_col=0)
size = jaccard.shape[0]

with open(out_path, 'w') as out_file:
    for i, node_a in enumerate(jaccard.index.values):
        for j, node_b in enumerate(jaccard.columns.values):
            if j >= i: continue
            out_file.write('\t'.join([str(node_a), str(node_b), str(jaccard.loc[node_a, node_b])]) + '\n')
