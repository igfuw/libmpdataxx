#!/usr/bin/env sh
set -e
cd tests/sandbox
mkdir build
cd build
cmake ..
# compiling and running pbl_iles on clang
# "/" intentional! (just to make cat exit with an error code)
if [[ $COMPILER == 'clang++' ]]; then VERBOSE=1 make pbl_iles_travis; fi
if [[ $COMPILER == 'clang++' ]]; then OMP_NUM_THREADS=4 ctest -V -R pbl_iles_travis || cat Testing/Temporary/LastTest.log /; fi
# smg on gcc
if [[ $COMPILER == 'g++' ]]; then VERBOSE=1 make pbl_smg_travis; fi
if [[ $COMPILER == 'g++' ]]; then OMP_NUM_THREADS=4 ctest -V -R pbl_smg_travis || cat Testing/Temporary/LastTest.log /; fi
cd ../../..
