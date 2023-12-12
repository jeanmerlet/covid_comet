import os, sys

# script arguments
input_path    = sys.argv[1]
output_dir    = sys.argv[2]
command_dir   = sys.argv[3]

# parameters
inflations = [1.2, 2.0, 4.0]
nodes = 100
wall_time = '6:00:00'

# make job script for each inflation value
for inflation in inflations:
    output_name = f'clusters_inflation-{inflation}.tsv'
    output_path = os.path.join(output_dir, output_name)
    command_name = f'jobscript_hipmcl_inf-{inflation}.sbatch'
    command_path = os.path.join(command_dir, command_name)

    lines = [
        '#!/bin/bash\n'
        f'#SBATCH -t {wall_time}\n',
        '#SBATCH -A SYB111\n',
        f'#SBATCH -N {nodes}\n',
        '#SBATCH -J hipmcl\n',
        f'#SBATCH -o ./logs/hipmcl_{inflation}.%J.out\n',
        f'#SBATCH -e ./logs/hipmcl_{inflation}.%J.err\n',
        '\n',
        'module load cray-mpich\n',
        '\n',
        'HIPMCL_EXE="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/hipmcl/bin/hipmcl"\n',
        'unset MPICH_OFI_NIC_POLICY\n',
        '\n',
        f'N={nodes}\n',
        'c=56\n',
        '\n',
        f'I={inflation}\n',
        'mem=450\n'
        '\n',
        f'IN_FILE="{input_path}"\n',
        f'OUT_FILE="{output_path}"\n',
        '\n',
        'srun -n $N -c $c $HIPMCL_EXE -M $IN_FILE -I $I -rand 0 -per-process-mem $mem -o $OUT_FILE\n'
    ]

    with open(command_path, 'w') as command_file:
        command_file.writelines(lines)
