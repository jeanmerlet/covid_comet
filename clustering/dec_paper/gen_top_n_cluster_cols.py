import pandas as pd
import numpy as np
import os, re

num_clusters = 'all'
clusters_filename = 'clusters_inflation-1.2.tsv'

#clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/frontier_run_out/threshold_3.1997/edge_list/hip_mcl'
clusters_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/thresh_run_HH-0.19991/hip_mcl'
clusters_path = os.path.join(clusters_dir, clusters_filename)
out_path = os.path.join(clusters_dir, 'top-' + str(num_clusters) + '_' + clusters_filename)

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

def get_top_country(meta, seq_ids):
    locations = meta.loc[seq_ids, 'Location'].values
    conts = []
    countries = []
    for location in locations:
        cont = location.split(' / ')[0]
        conts.append(cont)
        country = location.split(' / ')[1]
        countries.append(country)
    uniq_conts, uniq_conts_counts = np.unique(np.array(conts), return_counts=True)
    #top_cont = np.argsort(uniq_conts_counts)[::-1][0]
    top_cont = np.argsort(uniq_conts_counts)[::-1]
    conts = ','.join(uniq_conts[top_cont])
    conts_counts = ','.join([str(x) for x in uniq_conts_counts[top_cont]])
    uniq_countries, uniq_countries_counts = np.unique(np.array(countries), return_counts=True)
    top_country = np.argsort(uniq_countries_counts)[::-1][0]
    return conts, conts_counts, uniq_countries[top_country], uniq_countries_counts[top_country]
    #return uniq_conts[top_cont], uniq_conts_counts[top_cont], uniq_countries[top_country], uniq_countries_counts[top_country]

def get_sample_collection_date_median_and_range(meta, seq_ids):
    dates = meta.loc[seq_ids, 'Collection date'].values
    new_dates = []
    for date in dates:
        if date.split('-')[0] == '2019':
            new_dates.append('2020-01-01')
        elif len(date.split('-')) == 2:
            new_dates.append(date + '-' + '15')
        else:
            new_dates.append(date)
    new_dates = np.sort(new_dates)
    #if len(new_dates) == 0: return 0, '0', 0
    median_date = new_dates[round(len(new_dates) / 2)]
    med_year, med_month, med_day = median_date.split('-')
    median_date_num = ((int(med_year) - 2020) * 365) + (int(med_month) * 30) + int(med_day)
    earliest = new_dates[0].split('-')
    latest = new_dates[-1].split('-')
    earliest = [int(x) for x in earliest]
    latest = [int(x) for x in latest]
    time_diff = ((latest[0] - earliest[0]) * 365) + ((latest[1] - earliest[1]) * 30) + (latest[2] - earliest[2])
    return time_diff, median_date, median_date_num

# first find index of top 20 largest clusters
cluster_sizes = []
with open(clusters_path, 'rt') as clusters_file:
    for i, line in enumerate(clusters_file):
        seq_ids = line.strip().split(' ')
        num_seqs = len(seq_ids)
        cluster_sizes.append(num_seqs)

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
header = '\t'.join(['cluster_idx', 'cluster_size', 'date_range', 'median_date', 'median_days_from_01.01.20',
                    'top_cont', 'size_top_cont', 'top_country', 'size_top_country']) + '\n'
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
                if num_seqs > 1:
                    date_range, median_date, date_int = get_sample_collection_date_median_and_range(meta, seq_ids)
                    top_cont, size_top_cont, top_country, size_top_country = get_top_country(meta, seq_ids)
                    out_line = '\t'.join([str(i), str(num_seqs), str(date_range), median_date, str(date_int),
                                          top_cont, str(size_top_cont), top_country, str(size_top_country)]) + '\n'
                    out_file.write(out_line)
                    num_clusters_done += 1
                    if num_clusters != 'all' and num_clusters_done == num_clusters: break
            else:
                continue
