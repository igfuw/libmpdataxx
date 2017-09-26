#!/usr/bin/env sh
set -e
cd tests/unit
mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 $make_j
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /
cmake -DCMAKE_BUILD_TYPE=Release ../
VERBOSE=1 $make_j
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /
cd ../../..
