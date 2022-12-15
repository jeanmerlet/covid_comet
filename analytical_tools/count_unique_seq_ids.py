import os
import re
import numpy as np

root_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original'

fasta_names = []
for r, d, f in os.walk(root_dir):
    for fasta in f:
        if re.search('.*fasta$', fasta):
            if 'seq_from_each_10' not in fasta:
                fasta_names.append(os.path.join(r, fasta))

fasta_names.sort()

seq_sets = {}
for i, f in enumerate(fasta_names):
    print(i)
    seq_names = []
    with open(f) as in_file:
        for line in in_file:
            if '>' in line:
                seq_names.append(line)

    head, tail = os.path.split(f)
    seq_sets[tail] = seq_names

for i, item in enumerate(seq_sets.items()):
    print(i)
    seq_set, seq_names = item
    if i == 0:
        old_seq_names = seq_names
    else:
        old_seq_names = np.union1d(old_seq_names, seq_names)

num_uniq_seqs = len(old_seq_names)
print(num_uniq_seqs)
