import os
import re
import numpy as np
import pandas as pd
from mpi4py import MPI

seqs_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02'
out_dir = os.path.join(seqs_dir, 'uniq_ids')
meta_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/metadata/metadata_2022_06_02/metadata.tsv'

try:
    os.mkdir(out_dir)
except OSError:
    pass

fasta_names = []
for r, d, f in os.walk(seqs_dir):
    for fasta in f:
        if re.search('.*fasta$', fasta):
            fasta_names.append(os.path.join(r, fasta))

fasta_names.sort()

def get_non_uniq_headers(metadata):
    virus_names = metadata['Virus name']
    names, counts = np.unique(virus_names, return_counts=True)
    non_uniq_names = names[counts > 1]
    return non_uniq_names


#TODO: count how many aren't unique and save them somewhere
def fix_seq_headers(seq_set, out_path, non_uniq_names, acc_ids):
    non_uniq = False
    valid_headers = acc_ids.index.values
    with open(seq_set) as in_file:
        with open(out_path, 'w') as out_file:
            i = 0
            for line in in_file:
                if '>' in line:
                    i += 1
                    if i == 1:
                        out_file.write(line)
                        continue
                    name = line.split('|')[0][1:]
                    if name not in valid_headers:
                        non_uniq = True
                        continue
                    elif name in non_uniq_names:
                        non_uniq = True
                        continue
                    else:
                        non_uniq = False
                        name = '>' + '|'.join((name, acc_ids.loc[name, 'Accession ID'])) + '\n'
                        out_file.write(name)
                else:
                    if not non_uniq:
                        out_file.write(line)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, seq_set in enumerate(fasta_names):
    # distribute across ranks
    if i % size != rank: continue
    _, name = os.path.split(seq_set)
    out_path = os.path.join(out_dir, 'uniq_ids_' + name)
    metadata = pd.read_csv(meta_path, sep='\t', low_memory=False)
    non_uniq_names = get_non_uniq_headers(metadata)
    acc_ids = metadata[['Virus name', 'Accession ID']].set_index('Virus name')
    fix_seq_headers(seq_set, out_path, non_uniq_names, acc_ids)
