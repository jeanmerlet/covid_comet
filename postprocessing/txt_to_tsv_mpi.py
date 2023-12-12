from mpi4py import MPI
import argparse
import os, re


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
args = parser.parse_args()


in_file_list = [os.path.join(args.data_dir, f) for f in os.listdir(args.data_dir) if '.txt' in f]
in_file_list.sort()


def txt_to_tsv(in_path, out_path):
    with open(in_path, 'r') as in_file:
        for line in in_file:
            _, low_high1, _, low_high2, seq_id1, seq_id2, score = line.strip().split(' ')
            if low_high1 == '1' and low_high2 == '1':
                with open(out_path, 'a') as out_file:
                    out_file.write('\t'.join([seq_id1[:-2], seq_id2[:-2], score]) + '\n')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_file_path in enumerate(in_file_list):
    # match file to rank (1 per node)
    if i != rank: continue
    head, file_name = os.path.split(in_file_path)
    out_path = re.sub('.txt', '.tsv', in_file_path)
    print(f'Started {file_name} ({rank}/{size})...', flush=True)
    txt_to_tsv(in_file_path, out_path)
    print(f'Finished {file_name} ({rank}/{size}).', flush=True)
