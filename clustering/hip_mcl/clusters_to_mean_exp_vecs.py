import os, re, time
import pandas as pd
import numpy as np

clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990/edge_list/hip_mcl'
seq_data_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_01.tsv'

seq_ids_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/seq_ids.txt'
seq_ids = pd.read_csv(seq_ids_path, header=None, index_col=None, dtype=str).values
seq_row_idxs = {}
for i, seq_id in enumerate(seq_ids):
    seq_row_idxs[seq_id[0]] = i

clusters_paths = []
for r, d, f in os.walk(clusters_dir):
    for clusters_path in f:
        if 'variant' not in clusters_path:
            if 'vector' not in clusters_path:
                if 'intersections' not in clusters_path:
                    clusters_paths.append(os.path.join(r, clusters_path))
clusters_paths.sort()

def seq_ids_to_mean_vector(seq_ids, seq_data_path):
    idxs_to_read = []
    for i, seq_id in enumerate(seq_ids):
        idxs_to_read.append(seq_row_idxs[seq_id])
    cum_sum = np.zeros(73965)
    with open(seq_data_path) as seq_data:
        for i, line in enumerate(seq_data):
            if i not in idxs_to_read:
                continue
            else:
                seq = line.strip().split('\t')[4:]
                seq = np.array(seq, dtype=int)
                cum_sum = np.add(cum_sum, seq)

    start = time.time()
    num_seqs = len(seq_ids)
    mean_vector = cum_sum / num_seqs
    mean_vector = list(mean_vector)
    mean_vector = [str(round(x, 2)) for x in mean_vector]
    return mean_vector

total_time = time.time()
for i, clusters_path in enumerate(clusters_paths):
    _, name = os.path.split(clusters_path)
    inflation = re.search('\d\.\d', name).group(0)
    if inflation != '1.2': continue
    out_path = os.path.join(clusters_dir, 'mean_vectors-' + inflation)
    with open(clusters_path, 'r') as clusters:
        for j, cluster in enumerate(clusters):
            if j < 160: continue
            seq_ids = cluster.strip().split(' ')
            # hip_mcl output seq_ids have _A appended to the end
            seq_ids = [seq_id[:-2] for seq_id in seq_ids]
            single_cluster = time.time()
            mean_vector = seq_ids_to_mean_vector(seq_ids, seq_data_path)
            write_time = time.time()
            with open(out_path, 'a') as out_clusters:
                out_clusters.write('\t'.join(mean_vector) + '\n')
            print(f'inflation {inflation}: cluster {j} took {time.time() - single_cluster}', flush=True)

print(f'total time took: {time.time() - total_time}', flush=True)
