#!/bin/sh -e

set -x

QRY="${DATADIR}/probes.fwd.fa"
IN="${DATADIR}/probe2targets.expected.1.tsv"

opt -o $RESULTS stat -i $IN -q $QRY > "${RESULTS}/stat.out"

set +x

ACTUAL=$(wc -l < "${RESULTS}/stat_off_target_probes.txt")
TARGET=2
awk -v actual="$ACTUAL" -v target="$TARGET" ' 
    BEGIN {
        print (actual == target) ? "GOOD" : "BAD"
        print "Expected: " target;
        print "Actual: " actual;
    }' > "${RESULTS}.report"