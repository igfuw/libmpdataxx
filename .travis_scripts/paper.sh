#!/usr/bin/env sh
set -e
cd tests/paper_2015_GMD
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 make -C 6_coupled_harmosc
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make -C 6_coupled_harmosc test || cat 6_coupled_harmosc/Testing/Temporary/LastTest.log /

# compiling everything in the Release mode
cmake -DCMAKE_BUILD_TYPE=Release ../
VERBOSE=1 $make_j

# running all paper tests in Release mode 
#- OMP_NUM_THREADS=1 make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /
cd ../../..
