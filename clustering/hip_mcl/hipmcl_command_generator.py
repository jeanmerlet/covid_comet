import os, sys

# script arguments
input_path    = sys.argv[1]
output_dir    = sys.argv[2]
command_dir   = sys.argv[3]

# parameters
inflations = [1.2, 1.5, 2.0, 4.0, 8.0]
nodes = 25
wall_time = '2:00'

# make job script for each inflation value
for inflation in inflations:
    output_name = f'clusters_inflation-{inflation}.tsv'
    output_path = os.path.join(output_dir, output_name)
    command_name = f'submit_hipmcl_inflation-{inflation}.lsf'
    command_path = os.path.join(command_dir, command_name)

    lines = [
        '#!/bin/bash\n'
        f'#BSUB -W {wall_time}\n',
        '#BSUB -P SYB108\n',
        f'#BSUB -nnodes {nodes}\n',
        '#BSUB -J hipmcl\n',
        '#BSUB -o ./logs/hipmcl.%J.out\n',
        '#BSUB -e ./logs/hipmcl.%J.err\n',
        '\n',
        #'module load cuda/10.1.243 cmake/3.17.3 gcc/8.1.1 boost/1.70.0 spectrum-mpi/10.3.1.2-20200121\n',
        'module load cuda/10.1.243 cmake gcc boost\n',
        'module unload xalt\n',
        '\n',
        'HIPMCL_EXE=/gpfs/alpine/syb105/proj-shared/Personal/gabrielgaz/Apps/summit/hipMCL/gpu/bin/hipmcl-gpu\n',
        '\n',
        'export OMP_NUM_THREADS=4\n',
        '\n',
        '# -------------\n',
        '# JSRUN Configs\n',
        '# -------------\n',
        '\n',
        '# Total Resource Sets\n',
        f'N={nodes}\n',
        '# Resource Sets Per Node\n',
        'r=1\n',
        '# Number of Tasks - MPI Ranks\n',
        'a=1\n',
        '# Physical Cores per Resource Set\n',
        'c=42\n',
        '# GPUs Per Resource Set\n',
        'g=6\n',
        '\n',
        '# --------------\n',
        '# HipMCL Configs\n',
        '# --------------\n',
        '\n',
        '# Memory per Node - On Summit 512 GB / # of Processes Per Node\n',
        'mem=450\n',
        '# Inflation Value\n',
        f'I={inflation}\n',
        '\n',
        f'IN_FILE="{input_path}"\n',
        f'OUT_FILE="{output_path}"\n',
        '\n',
        'jsrun -n $N -r $r -a $a -c $c -g $g -bpacked:$c --smpiargs="-mca coll_ibm_skip_alltoallv true" $HIPMCL_EXE -M $IN_FILE -I $I -rand 0 -per-process-mem $mem -o $OUT_FILE'
    ]

    with open(command_path, 'w') as command_file:
        command_file.writelines(lines)
