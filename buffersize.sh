#!/usr/bin/env bash
NUM_RUNS=$1
CALL="${2}"
ARGS="${@:3}"
for ARG in $ARGS; do
    CALL="$CALL $ARG"
done
mkdir -p buffersize
TIMESTAMP=$(date +%d%m-%H%M%S)
echo "Buffer Size,Average Runtime" > buffersize/buffersize_$TIMESTAMP.csv
for i in {10..20}; do
    BUFFER_SIZE=$((2**$i))
    RUNTIME=0
    for RUN_NO in $(seq 1 $NUM_RUNS); do
        echo "RUN #$RUN_NO/$NUM_RUNS FOR BUFFER SIZE $BUFFER_SIZE" 
        START=$(date +%s)
        python $CALL -b $BUFFER_SIZE > /dev/null 
        END=$(date +%s)
        rm -r codes3d_summary
        ((RUNTIME += END-START))
    done
    AVG_RUN=$((RUNTIME/NUM_RUNS))
    echo "$BUFFER_SIZE,$AVG_RUN" >> buffersize/buffersize_$TIMESTAMP.csv
done
