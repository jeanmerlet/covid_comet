from mpi4py import MPI
import subprocess
import argparse
import re
import os

# data directory
data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/2-way-0.965_mut-100'
# tped filename prefix
tped_prefix = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/combined_na-250k'

# postprocessing tools binary dir
comet_tools_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/other_tools/comet_devel/genomics_gpu/tools'
num_way = '2'

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--start_idx')
parser.add_argument('-e', '--end_idx')
args = parser.parse_args()

def comet_postprocess(in_file):
    postprocess = os.path.join(comet_tools_dir, 'postprocess')
    allele_labels = tped_prefix + '_allele_labels.txt'
    line_labels = tped_prefix + '_line_labels.txt'
    out_file = re.sub(r'bin$', 'txt', in_file)
    
    subprocess.run([postprocess,
                    num_way,
                    allele_labels,
                    line_labels,
                    in_file,
                    out_file])

in_file_list = []
for r, d, f in os.walk(data_dir):
  for bin_file in f:
    if '.bin' in bin_file:
      in_file_list.append(os.path.join(r, bin_file))

in_file_list.sort()
in_file_list = in_file_list[int(args.start_idx):int(args.end_idx)]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_file_path in enumerate(in_file_list):
    # distribute across ranks
    if i % size != rank: continue
    head, file_name = os.path.split(in_file_path)
    print(f'postprocessing {file_name} ({rank}/{size})')
    comet_postprocess(in_file_path)
