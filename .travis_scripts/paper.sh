#!/usr/bin/env sh
set -e
cd tests/paper_2015_GMD
mkdir build
cd build
# Travis default is not the packaged one
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi
# the one from homebrew
if [[ $TRAVIS_OS_NAME == 'osx' && $CXX == 'g++' ]]; then cmake -DCMAKE_CXX_COMPILER=g++-4.8 ../; fi
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
if [[ $MPI == 'none' ]]; then OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /; fi

# some tests take too long with mpi, so we skip them
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C 0_basic_example test || cat 0_basic_example/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C 1_advscheme_opts test || cat 1_advscheme_opts/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C 6_coupled_harmosc test || cat 6_coupled_harmosc/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C 7_shallow_water test || cat 7_shallow_water/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C 8_boussinesq_2d test || cat 8_boussinesq_2d/Testing/Temporary/LastTest.log /; fi
cd ../../..
