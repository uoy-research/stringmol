#!/bin/bash



cd ../src
make clean
make debug

cd ../tests

#remove any previous working version
rm -f test

# Stop execution on any error
set -e

echo "building tests..."

RP="../debug"
#TODO: can't use wildcards e.g. ../release/*.o because of multiple 'main's... fix!
g++ -std=gnu++11 -Wall  -DDEBUG -g -o test  *.cpp ${RP}/mt19937-2.o ${RP}/randutil.o \
  ${RP}/agent.o ${RP}/SMspp.o ${RP}/stringPM.o ${RP}/rules.o \
  ${RP}/alignment.o ${RP}/params.o ${RP}/memoryutil.o ${RP}/stringmanip.o \
  ${RP}/hsort.o ${RP}/opcodes.o ${RP}/setupSM.o ${RP}/lodepng.o

echo "...success!"  
  
  
cd ../output
echo ""
echo "  now testing.."
../tests/test


echo "  cleaning up.."
rm -f rng.txt

cd ../tests


#############
# todo: consider using cccc to get software metrics..
# todo: consider how to use valgrind within this
