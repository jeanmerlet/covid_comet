import pandas as pd
import numpy as np
import os, re

num_clusters = 'all'
num_lineages = 'all'
clusters_filename = 'clusters_inflation-1.2.tsv'

#clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/frontier_run_out/threshold_3.1997/edge_list/hip_mcl'
clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/thresh_run_HH-0.19991/hip_mcl'
clusters_path = os.path.join(clusters_dir, clusters_filename)
out_path = os.path.join(clusters_dir, 'top-' + str(num_lineages) + '-lineages_in-top-' + str(num_clusters) + '-clusters_' + clusters_filename)

#metadata_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2021_09_30/metadata.tsv'
metadata_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/metadata.tsv'
meta = pd.read_csv(metadata_path, sep='\t', low_memory=False)
meta_seq_id_col = []
epi_isl_ids = meta.loc[:, 'Accession ID'].values
for seq_id in epi_isl_ids:
    new_id = seq_id.split('_')[-1].zfill(9)
    meta_seq_id_col.append(new_id)

meta.index = meta_seq_id_col

def format_ids(seq_ids):
    new_ids = []
    for seq_id in seq_ids:
        new_id = seq_id.split('_')[0]
        # omit wuhan refseq
        if new_id != '000000000':
            new_ids.append(new_id)
    return new_ids

def remove_year_only_ids(seq_ids):
    dates = meta.loc[seq_ids, 'Collection date'].values
    date_lens = []
    for date in dates:
        date_lens.append(len(date.split('-')))
    dates_with_month = np.array(date_lens) > 1
    dates_with_ids = meta.loc[seq_ids, 'Collection date']
    valid_ids = dates_with_ids.index.values[dates_with_month]
    return valid_ids

def get_top_n_lineages(n, meta, seq_ids):
    lineages = meta.loc[seq_ids, 'Pango lineage'].values
    lineages = lineages.astype('str')
    uniq, uniq_counts = np.unique(lineages, return_counts=True)
    if n != 'all':
        top_lineages = np.argsort(uniq_counts)[::-1][:n]
    else:
        top_lineages = np.argsort(uniq_counts)[::-1]
    top_n_lineages = uniq[top_lineages]
    top_n_lineages = ','.join(top_n_lineages)
    top_n_lineages_sizes = uniq_counts[top_lineages]
    top_n_lineages_sizes = ','.join([str(x) for x in top_n_lineages_sizes])
    return top_n_lineages, top_n_lineages_sizes

# sort by largest cluster
cluster_sizes = []
with open(clusters_path, 'rt') as clusters_file:
    for i, line in enumerate(clusters_file):
        seq_ids = line.strip().split(' ')
        num_seqs = len(seq_ids)
        cluster_sizes.append(num_seqs)

# subset to top n largest clusters
if num_clusters != 'all':
    top_n_idx = np.argsort(np.array(cluster_sizes))[::-1][:num_clusters]
else:
    top_n_idx = np.argsort(np.array(cluster_sizes))[::-1]

# get percent intersection of seq_ids with other clusters for each of the top n clusters
# this makes no sense within a particular inflation since the clusters are composed of unique seq_ids
seq_ids_by_cluster = {}
with open(clusters_path, 'rt') as clusters_file:
    for i, line in enumerate(clusters_file):
        break
        if i in top_n_idx:
            seq_ids = line.strip().split(' ')
            seq_ids = format_ids(seq_ids)
            seq_ids = remove_year_only_ids(seq_ids)
            seq_ids_by_cluster[i] = seq_ids

top_prop_intersect_id = {}
top_prop = {}
for cluster_id, seq_ids in seq_ids_by_cluster.items():
    break
    num_seqs = len(seq_ids)
    intersection_props, intersection_ids = [], []
    for other_cluster_id, other_seq_ids in seq_ids_by_cluster.items():
        if other_cluster_id == cluster_id:
            continue
        else:
            size_intersection = len(np.intersect1d(seq_ids, other_seq_ids))
            intersection_props.append(size_intersection / num_seqs)
            intersection_ids.append(other_cluster_id)
    top_id_idx = np.argsort(intersection_props)[-1]
    top_prop_intersect_id[cluster_id] = intersection_ids[top_id_idx]
    top_prop[cluster_id] = intersection_props[top_id_idx]

# get meta information for top n clusters
header = '\t'.join(['cluster_idx', 'cluster_size', 'top_lineages', 'size_top_lineages']) + '\n'
num_clusters_done = 0
cluster_ids = {}
with open(out_path, 'wt') as out_file:
    out_file.write(header)
    with open(clusters_path, 'rt') as clusters_file:
        for i, line in enumerate(clusters_file):
            if i in top_n_idx:
                seq_ids = line.strip().split(' ')
                seq_ids = format_ids(seq_ids)
                seq_ids = remove_year_only_ids(seq_ids)
                num_seqs = len(seq_ids)
                top_lineages, size_top_lineages = get_top_n_lineages(num_lineages, meta, seq_ids)
                out_line = '\t'.join([str(i), str(num_seqs), top_lineages, str(size_top_lineages)]) + '\n'
                out_file.write(out_line)
                num_clusters_done += 1
                if num_clusters != 'all' and num_clusters_done == num_clusters: break
            else:
                continue
