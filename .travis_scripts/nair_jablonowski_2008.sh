#!/usr/bin/env sh
set -e
cd tests/nair_jablonowski_2008
mkdir build
cd build
# the one from homebrew
if [[ $TRAVIS_OS_NAME == 'osx' && $CXX == 'g++' ]]; then cmake -DCMAKE_CXX_COMPILER=g++-4.8 ../; fi

# compiling everything in the Release mode
cmake -DCMAKE_BUILD_TYPE=Release ../
VERBOSE=1 $make_j

# running tests in Release mode 
# "/" intentional! (just to make cat exit with an error code)
OMP_NUM_THREADS=4 ctest -V -I 1,6 || cat Testing/Temporary/LastTest.log /
cd ../../..
