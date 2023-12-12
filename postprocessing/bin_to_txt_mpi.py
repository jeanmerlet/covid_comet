from mpi4py import MPI
import subprocess
import argparse
import re, os


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-t', '--tped_prefix')
parser.add_argument('-c', '--comet_tools_dir')
parser.add_argument('-w', '--num_way')
args = parser.parse_args()


def comet_postprocess(in_file, comet_tools_dir, tped_prefix, num_way):
    postprocess = os.path.join(comet_tools_dir, 'postprocess')
    allele_labels = tped_prefix + '_allele_labels.txt'
    line_labels = tped_prefix + '_line_labels.txt'
    out_file = re.sub(r'bin$', 'txt', in_file)
    subprocess.run([postprocess,
                    'duo',
                    num_way,
                    allele_labels,
                    line_labels,
                    in_file,
                    out_file])


in_file_list = [os.path.join(args.data_dir, f) for f in os.listdir(args.data_dir) if '.bin' in f]
in_file_list.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_file_path in enumerate(in_file_list):
    # distribute across ranks
    if i != rank: continue
    head, file_name = os.path.split(in_file_path)
    print(f'Postprocessing {file_name} ({rank}/{size})...', flush=True)
    comet_postprocess(in_file_path, args.comet_tools_dir, args.tped_prefix, args.num_way)
    print(f'Finished {file_name} ({rank}/{size}).', flush=True)
