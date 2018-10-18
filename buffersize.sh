#!/usr/bin/env bash
NUM_RUNS=$1
CALL="${2}"
ARGS="${@:3}"
for ARG in $ARGS; do
    CALL="$CALL $ARG"
done
for i in {10..20}; do
    RUNTIME=0
    for RUN_NO in $(seq 1 $NUM_RUNS); do
        BUFFER_SIZE=$((2**$i))
        START=$(date +%s)
        python $CALL -b $BUFFER_SIZE &> /dev/null #2> "error_logs/stderr_$RUN_NO"
        END=$(date +%s)
        rm -r codes3d_summary
        ((RUNTIME += END-START))
    done
    AVG_RUN=$((RUNTIME/NUM_RUNS))
    echo "BUFFER SIZE: $BUFFER_SIZE AVERAGE RUNTIME: ${AVG_RUN}s" 
done
