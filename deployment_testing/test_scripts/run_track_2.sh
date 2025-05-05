#!/bin/sh -e

set -x

QRY="${DATADIR}/probes.fwd.fa"
TGT_1="${DATADIR}/annotation.gff"
TGT_2="${DATADIR}/annotation.fa"

opt -o $RESULTS -p $THREADS track -q $QRY -a $TGT_1 -t $TGT_2 -1
"${SCRIPTS}/check_eq.py" "${DATADIR}/probe2targets.expected.2.tsv" \
    "${RESULTS}/probe2targets.tsv" "${RESULTS}/eq.out"

set +x

ACTUAL=$(cat "${RESULTS}/eq.out")
TARGET=1
awk -v actual="$ACTUAL" -v target="$TARGET" ' 
    BEGIN {
        print (actual == target) ? "GOOD" : "BAD"
        print "Expected: " target;
        print "Actual: " actual;
    }' > "${RESULTS}.report"