#!/usr/bin/env sh
set -e
cd tests/sandbox
mkdir build
cd build
# Travis default is not the packaged one
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi
cmake ..
# compiling and running pbl_iles on clang
# "/" intentional! (just to make cat exit with an error code)
if [[ $CXX == 'clang++' ]]; then VERBOSE=1 make pbl_iles_travis; fi
if [[ $CXX == 'clang++' ]]; then OMP_NUM_THREADS=4 ctest -V -R pbl_iles_travis || cat Testing/Temporary/LastTest.log /; fi
# smg on gcc
if [[ $CXX == 'g++' ]]; then VERBOSE=1 make pbl_smg_travis; fi
if [[ $CXX == 'g++' ]]; then OMP_NUM_THREADS=4 ctest -V -R pbl_smg_travis || cat Testing/Temporary/LastTest.log /; fi
cd ../../..
