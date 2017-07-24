#!/usr/bin/env sh
set -e
cd tests/unit
mkdir build 
cd build
# the one from homebrew
if [[ $TRAVIS_OS_NAME == 'osx' && $CXX == 'g++' ]]; then cmake -DCMAKE_CXX_COMPILER=g++-4.8 ../; fi
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 $make_j
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /
cd ../../..
