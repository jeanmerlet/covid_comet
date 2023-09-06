from mpi4py import MPI
import os, re

data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/3-way-0.999_mut-100'

paths = []
for r, d, f in os.walk(data_dir):
    for txt in f:
        if re.match('^out\_\d+\.txt$', txt):
            path = os.path.join(r, txt)
            paths.append(os.path.join(r, txt))

paths.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(f'size: {size}')

for i, path in enumerate(paths):
    # distribute ranks across nodes
    if i % size != rank: continue
    print(f'rank: {rank}')
    #start = time.time()
    #head, tail = os.path.split(path)
    #file_num = re.search('.*\_(\d+).txt', tail).groups()[0]
    #out_path = os.path.join(out_dir, 'deduped_' + file_num + '.tsv')
    #deduplicate(path, out_path, file_num)
    #print(f'deduplicated {tail} ({rank}/{size}). time: {time.time() - start}')
