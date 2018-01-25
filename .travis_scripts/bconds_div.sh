#!/usr/bin/env sh
set -e
cd tests/sandbox
mkdir build
cd build
cmake ..
VERBOSE=1 $make_j
# running bconds_div in Release mode
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make -C 6_bconds_div test || cat 6_bconds_div/Testing/Temporary/LastTest.log /
cd ../../..
