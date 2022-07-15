import os, re
import numpy as np

cluster_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl'
out_dir = os.path.join(cluster_dir, 'cyto_attributes')

all_paths = os.listdir(cluster_dir)
cluster_paths = []

for path in all_paths:
    if re.match('^clusters', path):
        cluster_paths.append(os.path.join(cluster_dir, path))
cluster_paths.sort()

num_clusters = 100
out_path = os.path.join(out_dir, f'top_{num_clusters}_cluster_sizes.tsv')
with open(out_path, 'w') as out_file:
    out_file.write('name\tcluster_size\n')

for path in cluster_paths:
    _, name = os.path.split(path)
    inf = re.search('\d\.\d', name).group(0)

    sizes = []
    all_cluster_sizes = []
    with open(path, 'r') as clusters_file:
        for i, line in enumerate(clusters_file):
            sizes.append(len(line.strip().split(' ')))
            all_cluster_sizes.append(str(len(line.strip().split(' '))))
    # get idx of largest n clusters
    top_n_idx = np.argsort(np.array(sizes))[::-1][:num_clusters]

    cluster_sizes = []
    cluster_names = []
    for i, idx in enumerate(top_n_idx):
        cluster_sizes.append(all_cluster_sizes[idx])
        cluster_names.append(str(i + 1) + '-' + inf)

    with open(out_path, 'a') as out_file:
        for i in range(num_clusters):
            out_file.write(cluster_names[i] + '\t' + cluster_sizes[i] + '\n')
