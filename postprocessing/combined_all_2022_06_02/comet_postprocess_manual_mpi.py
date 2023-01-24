from mpi4py import MPI
import subprocess
import os, re

# directory
data_dir = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/thresh_run_HH-0.19991'
# filename
tped_prefix = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/ns-as-0s_mutation_count_filtered_100/combined_all_2022_06_02_ns-0s'
comet_tools_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/other_tools/comet_devel/genomics_gpu/tools'
num_way = '2'

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


# load a txt list of full filenames here
names_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/scripts/postprocessing/combined_all_2022_06_02/jobs/manual.txt'
with open(names_path) as f:
    in_file_list = f.readlines()

in_file_list = [f.strip() for f in in_file_list]


# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_file_path in enumerate(in_file_list):
    # distribute across ranks
    if i % size != rank: continue
    head, file_name = os.path.split(in_file_path)
    print(f'postprocessing {file_name} ({rank}/{size})')
    comet_postprocess(in_file_path)
