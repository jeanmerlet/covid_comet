#!/bin/bash 
  
#SBATCH -A SYB105
#SBATCH -t 0:30:00
#SBATCH -J mut_counts
#SBATCH -o ./logs/mut_counts%J.out
#SBATCH -e ./logs/mut_counts%J.err
#SBATCH -N 13

module load python
source activate base

PYFILE="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/scripts/preprocessing/count_bp_mutations.py"

srun -n 52 -c 8 python $PYFILE
