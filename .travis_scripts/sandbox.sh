#!/usr/bin/env sh
set -e
cd tests/sandbox
mkdir build
cd build
cmake ..
VERBOSE=1 $make_j
# running selected sandbox tests in Release mode
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make -C convergence_2d_3d test || cat convergence_2d_3d/Testing/Temporary/LastTest.log /
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then OMP_NUM_THREADS=4 make -C convergence_spacetime test || cat convergence_spacetime/Testing/Temporary/LastTest.log /; fi
cd ../../..
