#!/bin/sh -e

set -x

QRY="${DATADIR}/probes.fa"
TGT_1="${DATADIR}/annotation.gff"
TGT_2="${DATADIR}/annotation.fa"

opt -o $RESULTS flip -q $QRY -a $TGT_1 -t $TGT_2 > "${RESULTS}/flip.out"

set +x

ACTUAL=$(wc -l < "${RESULTS}/rev_cmped_probes.txt")
TARGET=2
awk -v actual="$ACTUAL" -v target="$TARGET" ' 
    BEGIN {
        print (actual >= target) ? "GOOD" : "BAD"
        print "Expected: " target;
        print "Actual: " actual;
    }' > "${RESULTS}.report"