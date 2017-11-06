#!/usr/bin/env sh
set -e
cd tests/unit
mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 $make_j
# "/" intentional! (just to make cat exit with an error code)
if [[ $MPI == 'none' ]]; then OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /; fi

# don't run Debug bconds and absorber tests on MPI as it takes too long
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C kahan_sum test || cat kahan_sum/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C cone_bugs test || cat cone_bugs/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C shallow_water test || cat shallow_water/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C concurrent_1d test || cat concurrent_1d/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C git_revision test || cat git_revision/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C var_dt test || cat var_dt/Testing/Temporary/LastTest.log /; fi

cmake -DCMAKE_BUILD_TYPE=Release ../
# on MPI run absorber and bconds in release mode
if [[ $MPI != 'none' ]]; then VERBOSE=1 make -C absorber; fi
if [[ $MPI != 'none' ]]; then VERBOSE=1 make -C bconds; fi
if [[ $MPI != 'none' ]]; then make -C absorber test || cat absorber/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C bconds test || cat bconds/Testing/Temporary/LastTest.log /; fi

if [[ $MPI == 'none' ]]; then VERBOSE=1 $make_j; fi
# excluding test_issue because it (sometimes) fails in Release mode
# for unknown reasons and it only seems to happen on Travis ...
# "/" intentional! (just to make cat exit with an error code)
if [[ $MPI == 'none' ]]; then OMP_NUM_THREADS=4 ctest -E test_issue || cat Testing/Temporary/LastTest.log /;fi
