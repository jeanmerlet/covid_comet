from mpi4py import MPI
import numpy as np
import time
import os, re


data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100'
out_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100/dedup'


paths = []
for r, d, f in os.walk(data_dir):
    for txt in f:
        if re.match('^out\_[0-9\-]+\.txt$', txt):
            path = os.path.join(r, txt)
            paths.append(os.path.join(r, txt))

paths.sort()


def write_edges(out_path, edges):
    with open(out_path, 'wt') as out_file:
        for pair, weight in edges.items():
            out_file.write(pair + '\t' + weight + '\n')


def test_edges(edges, id1, id2, id3, weight, counter):
    id1, id2, id3 = np.sort([id1, id2, id3])
    pairs = [id1[:-2] + '\t' + id2[:-2],
             id1[:-2] + '\t' + id3[:-2],
             id2[:-2] + '\t' + id3[:-2]]
    for pair in pairs:
        if pair in edges:
            counter += 1
            if edges[pair] < weight:
                edges[pair] = weight
        else:
            edges[pair] = weight
    return edges, counter


def deduplicate(in_path, out_path, file_num, check_weight='max'):
    edges = {}
    counter = 0
    with open(in_path, 'rt') as in_file:
        for i, line in enumerate(in_file):
            _, _, _, _, _, _, id1, id2, id3, weight = line.strip().split(' ')
            edges, counter = test_edges(edges, id1, id2, id3, weight, counter)
    write_edges(out_path, edges)
    print(f'number of dups in {file_num}: {counter}')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, path in enumerate(paths):
    # distribute ranks across nodes
    if i % size != rank: continue
    head, tail = os.path.split(path)
    print(f'deduplicating {tail} ({rank + 1}/{size})...', flush=True)
    start = time.time()
    file_num = re.search('.*\_([0-9\-]+).txt', tail).groups()[0]
    out_path = os.path.join(out_dir, 'deduped_' + file_num + '.tsv')
    deduplicate(path, out_path, file_num)
    print(f'deduplicated {tail} ({rank + 1}/{size}). time: {time.time() - start}', flush=True)
