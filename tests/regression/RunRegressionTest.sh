#!/bin/bash
# Regression test for spatial stringmol (trial type 30).
#
# Runs the v2.0.3 spatial config with a fixed seed for 10,000 steps and checks the
# population dynamics (popdy001.dat) against a stored copy. Any change to the
# chemistry, the order of random-number draws or the order molecules are processed
# in shows up as a difference, so a refactor or performance change that passes this
# test has not changed the results.
#
# If a change is meant to alter the dynamics, regenerate the expected file with
#   ./RunRegressionTest.sh --update
# and say why in the commit message.
#
# Usage (from the tests directory, after building with `make` in src):
#   ./RunRegressionTest.sh [--update]

set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
BIN="$ROOT/release/stringmol"
EXPECTED="$HERE/expected_popdy_v2.0.3_seed4357_10000.dat"

[ -x "$BIN" ] || { echo "No binary at $BIN - run make in src first"; exit 2; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

sed -e "s/^RANDSEED.*/RANDSEED\t4357/" \
    -e "s/^REPORTEVERY.*/REPORTEVERY 1000/" \
    -e "s/^IMAGEEVERY.*/IMAGEEVERY\t100000000/" \
    -e "s/^NSTEPS.*/NSTEPS\t\t10001/" \
    "$ROOT/config/spatial_stringmol_v2.0.3.conf" > run.conf

echo "Running spatial stringmol for 10,000 steps (takes about a minute)..."
"$BIN" 30 run.conf > stdout.txt 2>&1

if [ "$1" == "--update" ]; then
    cp popdy001.dat "$EXPECTED"
    echo "Updated $EXPECTED"
    exit 0
fi

if cmp -s popdy001.dat "$EXPECTED"; then
    echo "PASS: population dynamics match $(basename "$EXPECTED")"
else
    first=$(paste -d'|' popdy001.dat "$EXPECTED" | awk -F'|' '$1!=$2{split($1,a,",");print a[1];exit}')
    echo "FAIL: population dynamics differ from $(basename "$EXPECTED") (first difference at step ${first:-?})"
    exit 1
fi
