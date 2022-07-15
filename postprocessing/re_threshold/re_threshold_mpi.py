from mpi4py import MPI
import subprocess
import re, os, time
import argparse

threshold = 3.1995

parser = argparse.ArgumentParser()
parser.add_argument('file_num_min', type=int)
parser.add_argument('file_num_max', type=int)
args = parser.parse_args()

# directory
data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/frontier_run_out'
out_dir = f'/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2021_09_30/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-1s_mutation_count_filtered_100/frontier_run_out/threshold_{threshold}'

all_paths = os.listdir(data_dir)
paths = []
for path in all_paths:
    if '.txt' in path:
        if int(re.search('_(\d+)\.txt', path).groups(1)[0]) >= args.file_num_min:
            if int(re.search('_(\d+)\.txt', path).groups(1)[0]) < args.file_num_max:
                paths.append(os.path.join(data_dir, path))

paths.sort()

def re_threshold(in_path, out_path, threshold):
    with open(out_path, 'wt') as out_file:
        with open(in_path, 'rt') as in_file:
            for line_num, line in enumerate(in_file):
                line_threshold = float(line.strip().split(' ')[-1])
                if line_threshold > threshold:
                    out_file.write(line)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_path in enumerate(paths):
    if i == 0: start = time.time()
    # distribute across ranks
    if i % size != rank: continue
    head, name = os.path.split(in_path)
    print(f'postprocessing {name} ({rank}/{size})')
    out_path = os.path.join(out_dir, name)
    re_threshold(in_path, out_path, threshold)
    if i == 0: print(f'total time: {time.time() - start}')
