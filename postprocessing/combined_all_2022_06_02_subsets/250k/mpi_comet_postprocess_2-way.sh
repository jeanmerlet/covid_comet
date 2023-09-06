#!/bin/bash

COMET_TOOLS_DIR="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/other_tools/comet_devel/genomics_gpu/tools"
export PATH="${PATH}:$COMET_TOOLS_DIR"

# uncomment these 3 lines if recompilation is needed (andes -> summit for example)
pushd $COMET_TOOLS_DIR
./make.sh
popd

export PY_FILE="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sars_cov_2_geo/scripts/postprocessing/combined_all_2022_06_02_subsets/250k/comet_postprocess_single_mpi_250k_2-way.py"
export COMET_OUT_FILES_DIR="/gpfs/alpine/syb105/proj-shared/Projects/GeoBio_CoMet/data/aligned/sequences_2022_06_02/uniq_ids/preprocessed_d-cutoff_1000_n-cutoff_0.01_pos_342-29665/usa_mex_can_250k/ns-as-0s_mutation_count_filtered_100/2-way-no-thresh_mut-100"

export NUM_FILES=$(find $COMET_OUT_FILES_DIR -name '2_way_*.bin' |  tr ' ' '\n' | wc -l)
export NUM_NODES_PER_BATCH=50
# correct rounding for truncating arithmetic
export NUM_BATCHES=$((( $NUM_FILES + $NUM_NODES_PER_BATCH - 1) / $NUM_NODES_PER_BATCH ))

for BATCH in $(seq 1 $NUM_BATCHES) ; do

    export CURRENT_BATCH=$BATCH
    bsub \
        -P syb111 \
        -nnodes $NUM_NODES_PER_BATCH \
        -W 6:00 \
        -J comet_post \
        -o ./logs/comet_post_2-way.%J.out \
        -e ./logs/comet_post_2-way.%J.err \
        ./comet_srun.sh

done
