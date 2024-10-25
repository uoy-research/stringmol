#!/bin/bash


#g++ main.cpp -o travis_gcc_cpp98
#./travis_gcc_cpp98
#cppcheck --quiet --error-exitcode=1 main.cpp

echo "============================"
echo "Checking with cppcheck"
cd src
# the --check-config argument seems to find standard library headers
#cppcheck --error-exitcode=1 --force --check-config --enable=all *.cpp
echo "----------------------------"
echo "Checking config with cppcheck"
cppcheck --error-exitcode=1 --force --check-config --suppress=missingIncludeSystem -i lodepng.* .
echo "----------------------------"
echo "Checking code with cppcheck"
cppcheck --error-exitcode=1 --force --enable=all --inline-suppr --suppress=missingIncludeSystem -i lodepng.* .
cd ../
echo ""

# NB: use "// cppcheck-suppress unusedFunction" before functions to suppress warnings
#           //cppcheck-suppress invalidscanf_libc"


echo "======================================"
echo "Running Smoke Tests"
echo "---------------------------------"
echo "Building!"


# Stop execution on any error
set -e

cd src
make clean
make all
cd ../output


echo "---------------------------------"
echo "Checking TTYPE 0 (1 on 1) runs ok"
../release/stringmol 0 ../config/quick_test0.conf > tmp_stdout.txt
echo "exit status is $?"
rm tmp_stdout.txt
sh ../util/rm_runfiles.sh
echo "---------------------------------"
echo "Checking TTYPE 1 (ALIFE XII) runs ok"
echo "exit status is $?"
../release/stringmol 1 ../config/quick_test1.conf > tmp_stdout.txt
rm tmp_stdout.txt
sh ../util/rm_runfiles.sh
echo "---------------------------------"
echo "Checking TTYPE 33 (Spatial Stringmol) runs ok"
echo "exit status is $?"
../release/stringmol 33 ../config/quick_test33.conf > tmp_stdout.txt
rm tmp_stdout.txt
sh ../util/rm_runfiles.sh
cd ../

# unStop execution on any error
#unset -e TODO: this doesn't recognise the "-e" option

GREEN='\033[0;32m'
NC='\033[0m' # No Color
echo "${GREEN}===============================================================================${NC}"
echo "Running Catch.hpp Tests.  Please Wait."
#g++ -Wall Shapes-Catch-Testing-Example/Source/Shapes-Catch-Testing-Example.cpp Shapes-Catch-Testing-Example/Source/Implementation/*.cpp -o main
#g++ -Wall Shapes-Catch-Testing-Example/Test/*.cpp Shapes-Catch-Testing-Example/Source/Implementation/*.cpp -o test

echo "  compiling..."

cd src
make clean
make debug

cd ../tests
RP="../debug"
#TODO: can't use wildcards e.g. ../release/*.o because of multiple 'main's... fix!
g++ -std=gnu++11 -Wall  -DDEBUG -g -o test  *.cpp ${RP}/mt19937-2.o ${RP}/randutil.o \
  ${RP}/SMspp.o ${RP}/stringPM.o ${RP}/agents_base.o ${RP}/rules.o ${RP}/alignment.o \
  ${RP}/params.o ${RP}/memoryutil.o ${RP}/instructions.o ${RP}/stringmanip.o \
  ${RP}/hsort.o
cd ..
echo ""
echo "  now testing.."
./tests/test


echo "  cleaning up.."
rm -f rng.txt
