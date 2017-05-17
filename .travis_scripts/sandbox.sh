#!/usr/bin/env sh
set -e
cd tests/sandbox
mkdir build
cd build
# Travis default is not the packaged one
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi
# the one from homebrew
if [[ $TRAVIS_OS_NAME == 'osx' && $CXX == 'g++' ]]; then cmake -DCMAKE_CXX_COMPILER=g++-4.8 ../; fi
cmake ..
VERBOSE=1 $make_j
# running selected sandbox tests in Release mode
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make -C 0_mpi_adv test || cat 0_mpi_adv/Testing/Temporary/LastTest.log /
OMP_NUM_THREADS=4 make -C 3_convergence_2d_3d test || cat 3_convergence_2d_3d/Testing/Temporary/LastTest.log /
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then OMP_NUM_THREADS=4 make -C 5_convergence_spacetime test || cat 5_convergence_spacetime/Testing/Temporary/LastTest.log /; fi
cd ../../..
