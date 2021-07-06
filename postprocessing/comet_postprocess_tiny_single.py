import subprocess
import argparse
import re
import os

in_file_path = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/mutation_count_filtered_100/small_combined/comet_output_seq_from_each_2/out_0.bin'
tped_prefix = '/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/original/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/mutation_count_filtered_100/small_combined/seq_from_each_2'
comet_tools_dir = "/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/other_tools/comet/genomics_gpu/tools"
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

comet_postprocess(in_file_path)
