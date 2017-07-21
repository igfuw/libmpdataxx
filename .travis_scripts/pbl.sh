#!/usr/bin/bash
set -e
cd tests/sandbox
mkdir build
cd build
# Travis default is not the packaged one
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi
cmake ..
# compiling and running pbl_iles on clang
# "/" intentional! (just to make cat exit with an error code)
if [[ $CXX == 'clang++' ]]; then VERBOSE=1 make pbl_iles; fi
if [[ $CXX == 'clang++' ]]; then OMP_NUM_THREADS=4 ctest -V -R pbl_iles || cat /Testing/Temporary/LastTest.log /; fi
# smg on gcc
if [[ $CXX == 'g++' ]]; then VERBOSE=1 make pbl_smg; fi
if [[ $CXX == 'g++' ]]; then OMP_NUM_THREADS=4 ctest -V -R pbl_smg || cat /Testing/Temporary/LastTest.log /; fi
cd ../../..
