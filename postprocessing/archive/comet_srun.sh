#!/bin/bash

start_idx=$(( $NUM_NODES_PER_BATCH * ($CURRENT_BATCH - 1) ))
end_idx=$(( $NUM_NODES_PER_BATCH * $CURRENT_BATCH ))

srun -N $NUM_NODES_PER_BATCH -n $NUM_NODES_PER_BATCH --ntasks-per-node=1 -c 8 python $PY_FILE -s $start_idx -e $end_idx
