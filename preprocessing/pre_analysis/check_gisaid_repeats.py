import os
import re
import numpy as np

root_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data'

fasta_names = []
for r, d, f in os.walk(root_dir):
    for fasta in f:
        if re.search('.*fasta$', fasta):
            fasta_names.append(os.path.join(r, fasta))

fasta_names.sort()

seq_sets = {}
for f in fasta_names:
    seq_names = []
    with open(f) as in_file:
        for line in in_file:
            if '>' in line:
                seq_names.append(line)

    head, tail = os.path.split(f)
    seq_sets[tail] = seq_names

#print(list(seq_sets.keys()))
intersection_names = []
for seq_set1, seq_names1 in seq_sets.items():
    for seq_set2, seq_names2 in seq_sets.items():
        if seq_set1 != seq_set2:
            num_intersections = len(np.intersect1d(seq_names1, seq_names2))
            if num_intersections > 1:
                print(f'name of set 1: {seq_set1} ({len(seq_names1)} seqs)')
                print(f'name of set 2: {seq_set2} ({len(seq_names2)} seqs)')
                print(f'# of repeats: {num_intersections}\n')
            elif num_intersections == 0:
                print(f'name of set 1: {seq_set1} ({len(seq_names1)} seqs)')
                print(f'name of set 2: {seq_set2} ({len(seq_names2)} seqs)')
                print(f'NO INTERSECTION\n')
                

