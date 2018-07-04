#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
BASE_DIR=$SCRIPT_DIR
DMFT_EXEC=$BASE_DIR/DMFT.py
TEST_DIR=$BASE_DIR/tests

allcount="BLAH"
failcount=0

cd $TEST_DIR

for file in $TEST_DIR/*.in; do
    echo Found test suite $file
    tests=$(grep "^#\!TEST" "$file")
    while read line; do
	let allcount++
	toks=($line)
	testid=${toks[1]}
	testname=$(sed "s/_/ /g" <<< $testid)
	printf "Running test: %-40s" "$testname..."
	runline="$DMFT_EXEC $file General.FileNamePrefix=$testid ${toks[@]:2}"
	if $runline 2>/dev/null; then
	    echo -e "[\e[00;32m PASSED \e[00m]"
	else
	    echo -e "[\e[00;31m FAILED \e[00m]"
	    let failcount++
	    $runline
	fi
    done <<<"$tests"
done

rm -f _*.hdf5

if [ $failcount -eq 0 ]; then
    echo -e "Complete: \e[00;32mPassed $allcount tests.\e[00m"
else
    echo -e "Complete: \e[00;31mFailed $failcount of $allcount tests.\e[00m "
fi
