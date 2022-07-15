import os
import pandas as pd

cluster_identity_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl'

paths = os.listdir(cluster_identity_dir)

ident_paths = []
for path in paths:
    if 'all_variant' in path:
        ident_paths.append(os.path.join(cluster_identity_dir, path))
ident_paths.sort()

num_variants = 3
num_clusters = 20
for path in ident_paths:
    top_n_variants = pd.read_csv(path, sep='\t', usecols=list(range(num_variants)), nrows=num_clusters, index_col=None, header=None)
    print(top_n_variants)
    break

