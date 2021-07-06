from mpi4py import MPI
import subprocess
import argparse
import re
import os

data_dir = "/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/count_filtered_0.01N_1000D/trimmed_ends_first_primer_overlap/mutation_count_filtered_100/comet_tped/combined/comet_output"
tped_prefix = "/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/count_filtered_0.01N_1000D/trimmed_ends_first_primer_overlap/mutation_count_filtered_100/comet_tped/combined/combined_all"
comet_tools_dir = "/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/other_tools/comet/genomics_gpu/tools"
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