#!/bin/bash

start_idx=$(( $NUM_NODES_PER_BATCH * ($CURRENT_BATCH - 1) ))
end_idx=$(( $NUM_NODES_PER_BATCH * $CURRENT_BATCH ))

jsrun -n $NUM_NODES_PER_BATCH -a 1 -c 20 python $PY_FILE -s $start_idx -e $end_idx
