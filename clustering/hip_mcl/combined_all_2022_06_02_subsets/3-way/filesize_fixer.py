import subprocess
import shlex
import os, re

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100'

paths = []
for r, d, f in os.walk(data_dir):
    for txt in f:
        if re.match('^out\_\d+\.txt$', txt):
            path = os.path.join(r, txt)
            paths.append(os.path.join(r, txt))

paths.sort()


def split_file(path):
    head, tail = os.path.split(path)
    prefix = os.path.join(head, tail[:-4] + '-')
    subprocess.run(['split', '-C 10000m', '-d', '--additional-suffix=.txt', path, prefix])


def remove_file(path):
    subprocess.run(['rm', path])


# 10gb to bytes
max_filesize = 10000000000

for path in paths:
    if os.path.getsize(path) > max_filesize:
        head, tail = os.path.split(path)
        print(f'splitting {tail}...', end='', flush=True)
        split_file(path)
        remove_file(path)
        print('done!', flush=True)
