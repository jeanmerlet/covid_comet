import time
import os, re

tsv_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100/dedup/combined_dedup_3-way.tsv'
out_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100/dedup/rededuped_combined_dedup_3-way.tsv'


def write_edges(out_path, edges):
    with open(out_path, 'wt') as out_file:
        for pair, weight in edges.items():
            out_file.write(pair + '\t' + weight + '\n')


def test_edges(edges, id1, id2, weight, counter):
    pair = id1 + '\t' + id2
    if pair in edges:
        counter += 1
        if edges[pair] < weight:
            edges[pair] = weight
    else:
        edges[pair] = weight
    return edges, counter


def deduplicate(in_path, out_path):
    edges = {}
    counter = 0
    with open(in_path, 'rt') as in_file:
        for i, line in enumerate(in_file):
            id1, id2, weight = line.strip().split('\t')
            edges, counter = test_edges(edges, id1, id2, weight, counter)
    write_edges(out_path, edges)
    print(f'number of dups in {file_num}: {counter}')


print('\n\nStarting deduplication of combined file...')
start = time.time()
deduplicate(tsv_path, out_path)
print(f'total time based on {file_num}: {time.time() - start}')
