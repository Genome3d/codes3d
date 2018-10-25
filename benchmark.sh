#!/usr/bin/env bash
if [ "${1}" = "--temp" ]
    then 
        TEMP=true
        NUM_RUNS=$2
        FILENAME=$(basename $3)
        CALL="${3}"
        ARGS="${@:4}"
    else 
        TEMP=false
        NUM_RUNS=$1
        FILENAME=$(basename $2)
        CALL="${2}"
        ARGS="${@:3}"
fi 
RUNTIME=0
for ARG in $ARGS; do
    CALL="$CALL $ARG"
done
TEMP_FILES=""
#mkdir error_logs
for RUN_NO in $(seq 1 $NUM_RUNS); do
    echo "BENCHMARK RUN #$RUN_NO/$NUM_RUNS"
    TEMP_FILES="$TEMP_FILES temp_$RUN_NO.log"
    START=$(date +%s)
    python $CALL > /dev/null
    END=$(date +%s)
    ((RUNTIME += END-START))
    rm -r codes3d_summary

    python benchmark.py --preprocess codes3d/codes3d.py
    python /usr/local/lib/python2.7/dist-packages/kernprof.py -l $CALL &> /dev/null 
    rm -r codes3d_summary
    python /usr/local/lib/python2.7/dist-packages/line_profiler.py "${FILENAME}.lprof" > "temp_$RUN_NO.log";
    rm "${FILENAME%.py}.py.lprof"
    python benchmark.py --postprocess codes3d/codes3d.py

done
if [ $TEMP = "false" ] 
    then 
        echo "PROCESSING BENCHMARK RESULTS"
        python benchmark.py --preprocess codes3d/codes3d.py
        python benchmark.py --compile $TEMP_FILES 
        python benchmark.py --postprocess codes3d/codes3d.py
        echo "PROCESSED BENCHMARK RESULTS"
    else
        :
fi
rm $TEMP_FILES
AVG_RUN=$((RUNTIME/NUM_RUNS))
echo "AVERAGE RUNTIME: ${AVG_RUN}s" 
