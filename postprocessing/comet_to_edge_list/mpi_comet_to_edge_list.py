import subprocess
import os
from mpi4py import MPI

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/combined_all_LL-3.1990'
out_dir = os.path.join(data_dir, 'edge_list')

try:
    os.mkdir(out_dir)
except OSError:
    pass

data_paths = []
for r, d, f in os.walk(data_dir):
    for file_path in f:
        if '.txt' in file_path:
            data_paths.append(os.path.join(r, file_path))

def convert_to_edge_list_tsv(path, out_dir):
    _, file_name = os.path.split(path)
    out_path = os.path.join(out_dir, file_name[:-4] + '.tsv')
    subprocess.call(f"tr ' ' '\t' < {path} | cut -f 5-7 > {out_path}", shell=True)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, path in enumerate(data_paths):
    # distribute across ranks
    if i % size != rank: continue
    convert_to_edge_list_tsv(path, out_dir)
