#!/bin/bash

#COMET_TOOLS_DIR="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/other_tools/comet/genomics_gpu/tools"
#pushd $COMET_TOOLS_DIR
#./make.sh
#popd

export PATH="${PATH}:$COMET_TOOLS_DIR"

export PY_FILE="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/scripts/postprocessing/comet_postprocess_single_mpi.py"
export COMET_OUT_FILES_DIR="/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/comet_tped/combined/comet_output"

export NUM_FILES=$(find $COMET_OUT_FILES_DIR -name 'out*.bin' |  tr ' ' '\n' | wc -l)
export NUM_NODES_PER_BATCH=20
export NUM_BATCHES=$(( $NUM_FILES / $NUM_NODES_PER_BATCH ))

for BATCH in $(seq 1 $NUM_BATCHES) ; do

    export CURRENT_BATCH=$BATCH
    sbatch \
        -A syb105 \
        -N $NUM_NODES_PER_BATCH \
        -t 36:00:00 \
        -J comet_post \
        -o ./logs/comet_post_single.o%J \
        -e ./logs/comet_post_single.e%J \
        ./comet_srun.sh

done
