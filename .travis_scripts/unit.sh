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
# excluding test_issue because it (sometimes) fails in Release mode
# for unknown reasons and it only seems to happen on Travis ...
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 ctest -E test_issue || cat Testing/Temporary/LastTest.log /
cd ../../..
