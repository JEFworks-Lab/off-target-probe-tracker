#!/bin/sh

SCRATCH=$1
THREADS=$2
mkdir -p $SCRATCH

export SCRIPTS="./test_scripts"
export DATADIR="./test_data"
export THREADS

# unzip data files
echo "unzipping test data..."
gzip -d "$DATADIR"/*.gz

TESTS=""
run_test() {
    if [ "$#" -lt 2 ]; then
        echo "A test requires >= 2 params! exiting..."
        exit 1
    fi
    NAME="$1"
    FILE="$2"
    shift
    shift
    TESTS="${TESTS} ${NAME}"
    export RESULTS="${SCRATCH}/${NAME}"
    mkdir -p $RESULTS
    START="$(date +%s)"
    "${SCRIPTS}/${FILE}" "$@"
    STATUS="$?"
    END="$(date +%s)"
    if [ "${STATUS}" = "0" ]; then
        if [ -f "${RESULTS}.report" ] && [ "$(echo $(head -n 1 "${RESULTS}.report"))" = "GOOD" ]; then
            rm -rf "${RESULTS}"
            # echo "blah"
        fi
    fi
    eval "${NAME}_TIME"="$((END-START))"
}

set +e

run_test FLIP "run_flip.sh"
run_test TRACK_1 "run_track_1.sh"
run_test TRACK_2 "run_track_2.sh"
run_test TRACK_3 "run_track_3.sh"
run_test STAT "run_stat.sh"

set -e
printf "\n"
ERR=0
for i in ${TESTS}; do
    VAL="${i}_TIME"
    eval TIME="\$$VAL"
    printf "\033[1m$i (Time: %ss)\033[0m\n" "${TIME}"
    STATUS="$(head -n 1 "${SCRATCH}/${i}.report")"
    if [ "$STATUS" != "GOOD" ]; then
        printf "\033[31mTEST FAILED\033[0m\n"
        ERR=$((ERR+1))
    else
        printf "\033[32mTEST SUCCESS\033[0m\n"
    fi
    cat "${SCRATCH}/${i}.report"
    printf "\n"
done

echo "cleaning up..."
gzip "$DATADIR"/*.fa
gzip "$DATADIR"/*.gff
rm "$DATADIR"/*.fxi