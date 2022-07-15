from mpi4py import MPI
import os, re

data_dir = '/gpfs/alpine/syb105/proj-shared/Personal/wellsag/comet_2way_output'
out_dir = os.path.join(data_dir, 'edge_list')

try:
    os.mkdir(out_dir)
except OSError:
    pass

out1, out2 = 'T', 'T'
out_type = 'HH'

data_paths = [os.path.join(data_dir, path) for path in os.listdir(data_dir) if '.txt' in path]
data_paths.sort()

def check_correct_duo_type(id1, id2):
    id1_type = re.search('.*_([AT])$', id1).groups(0)[0]
    id2_type = re.search('.*_([AT])$', id2).groups(0)[0]
    if out1 == id1_type and out2 == id2_type:
        return True
    else:
        return False
    

def convert_to_edge_list_tsv(in_path):
    _, file_name = os.path.split(path)
    out_path = os.path.join(out_dir, out_type + '_' + file_name[:-4] + '.tsv')
    with open(out_path, 'wt') as out_file:
        with open(in_path, 'rt') as in_file:
            for line in in_file:
                _, _, _, _, id1, id2, edge_weight = line.strip('\n').split(' ')
                if check_correct_duo_type(id1, id2):
                    out_file.write('\t'.join((id1, id2, edge_weight)) + '\n')
                    
                    
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, path in enumerate(data_paths):
    # distribute across ranks
    if i % size != rank: continue
    convert_to_edge_list_tsv(path)
