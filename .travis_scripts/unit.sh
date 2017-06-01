#!/usr/bin/env sh
set -e
cd tests/unit
mkdir build 
cd build
# Travis default is not the packaged one
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi
# the one from homebrew
if [[ $TRAVIS_OS_NAME == 'osx' && $CXX == 'g++' ]]; then cmake -DCMAKE_CXX_COMPILER=g++-4.8 ../; fi
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 $make_j
# "/" intentional! (just to make cat exit with an error code)
if [[ $MPI == 'none' ]]; then OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /; fi

# don't run Debug bconds test on MPI as it takes too long
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C kahan_sum test || cat kahan_sum/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C cone_bugs test || cat cone_bugs/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C shallow_water test || cat shallow_water/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C concurrent_1d test || cat concurrent_1d/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C git_revision test || cat git_revision/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C absorber test || cat absorber/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C var_dt test || cat var_dt/Testing/Temporary/LastTest.log /; fi
# on MPI run bconds in release mode
cmake -DCMAKE_BUILD_TYPE=Release ../
if [[ $MPI != 'none' ]]; then VERBOSE=1 make -C bconds /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C bconds test || cat bconds/Testing/Temporary/LastTest.log /; fi
cd ../../..
