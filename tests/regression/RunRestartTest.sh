#!/bin/bash
# Restart test for spatial stringmol (trial type 30) - issue #46.
#
# Runs the v2.0.3 spatial config with a fixed seed to step 7,000, then restarts a
# second run from the step-5,000 checkpoint (out1_05000.conf + RNGstate_5000.txt)
# and checks that both runs have the same population counts at steps 6,000 and
# 7,000. Species ids are not compared: a species that goes extinct before the
# checkpoint and later re-appears is given a new id in the restarted run, because
# checkpoints only record extant species.
#
# Usage (from the tests directory, after building with `make` in src):
#   ./RunRestartTest.sh

set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
BIN="$ROOT/release/stringmol"

[ -x "$BIN" ] || { echo "No binary at $BIN - run make in src first"; exit 2; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/full" "$WORK/resume"

cd "$WORK/full"
sed -e "s/^RANDSEED.*/RANDSEED\t4357/" \
    -e "s/^REPORTEVERY.*/REPORTEVERY 1000/" \
    -e "s/^IMAGEEVERY.*/IMAGEEVERY\t100000000/" \
    -e "s/^NSTEPS.*/NSTEPS\t\t7001/" \
    "$ROOT/config/spatial_stringmol_v2.0.3.conf" > run.conf
echo "Running to step 7,000..."
"$BIN" 30 run.conf > stdout.txt 2>&1

cp out1_05000.conf "$WORK/resume/run.conf"
cp RNGstate_5000.txt "$WORK/resume/"
cd "$WORK/resume"
echo "Restarting from the step 5,000 checkpoint..."
"$BIN" 30 run.conf > stdout.txt 2>&1

# per-step sorted species counts, for steps after the checkpoint
counts() { awk -F, '$1>5000{print $1, $3}' "$1" | sort -k1,1n -k2,2n; }
if diff <(counts "$WORK/full/popdy001.dat") <(counts "$WORK/resume/popdy001.dat") >/dev/null; then
    echo "PASS: restarted run reproduces the population counts at steps 6,000 and 7,000"
else
    echo "FAIL: restarted run differs from the original after the checkpoint"
    exit 1
fi
