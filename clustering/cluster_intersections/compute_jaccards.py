import os, re, time
import pandas as pd
import numpy as np

clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl'
out_dir = os.path.join(clusters_dir, 'jaccard')

file_names = os.listdir(clusters_dir)
clusters_paths = [os.path.join(clusters_dir, name) for name in file_names if re.match('^clusters', name)]
clusters_paths.sort()

num_inflations = len(clusters_paths)
num_cluster_comparisons = 100

cluster_sets = []
inflation_names = []
for inf, clusters_path in enumerate(clusters_paths):
    _, name = os.path.split(clusters_path)
    inflation = re.search('\d\.\d', name).group(0)
    inflation_names.append(inflation)
    cluster_sets.append([])
    # clusters are not ordered by largest size
    # first pass-through to find idx of largest n
    top_n_idx = []
    clusters = []
    sizes = []
    with open(clusters_path, 'r') as clusters_file:
        for i, line in enumerate(clusters_file):
            clusters.append(line.strip().split(' '))
            sizes.append(len(line.strip().split(' ')))
    top_n_idx = np.argsort(np.array(sizes))[::-1][:num_cluster_comparisons]

    for idx in top_n_idx:
        cluster_sets[inf].append(clusters[idx])

# all comparisons
names = []
for i in range(num_cluster_comparisons):
    for inf in inflation_names:
        names.append(f'{i+1}-{inf}')

size = num_inflations * num_cluster_comparisons
comparisons = pd.DataFrame(np.zeros((size, size)), index=names, columns=names)

# row block
for i in range(num_cluster_comparisons):
    print(i, flush=True)
    # row increment
    for n, cluster_set_a in enumerate(cluster_sets):
        start_cluster = cluster_set_a[i]
        # column block
        for j in range(num_cluster_comparisons):
            # column increment
            for m, cluster_set_b in enumerate(cluster_sets):
                target_cluster = cluster_set_b[j]
                if len(target_cluster) > len(start_cluster):
                    num_intersection = np.count_nonzero(np.in1d(target_cluster, start_cluster))
                else:
                    num_intersection = np.count_nonzero(np.in1d(start_cluster, target_cluster))
                num_union = np.union1d(start_cluster, target_cluster).size
                jaccard = num_intersection / num_union
                comparisons.iloc[i*num_inflations + n, j*num_inflations + m] = round(jaccard, 2)

out_path = os.path.join(out_dir, f'intersections_top_{num_cluster_comparisons}_clusters.tsv')
comparisons.to_csv(out_path, sep='\t', index=True, header=True)
