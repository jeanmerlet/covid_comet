#!/bin/bash 
  
#SBATCH -A SYB105
#SBATCH -t 1:00:00
#SBATCH -J mut_counts
#SBATCH -o ./logs/filter_by_mut%J.out
#SBATCH -e ./logs/filter_by_mut%J.err
#SBATCH -N 13

module load python
source activate base

PYFILE="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/scripts/preprocessing/filter_by_bp_mutation.py"

srun -n 52 -c 8 python $PYFILE
