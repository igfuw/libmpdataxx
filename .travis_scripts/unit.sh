#!/usr/bin/env sh
set -e
cd tests/unit
mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../
VERBOSE=1 $make_j
# "/" intentional! (just to make cat exit with an error code)
if [[ $MPI == 'none' ]]; then OMP_NUM_THREADS=4 make test || cat Testing/Temporary/LastTest.log /; fi

# don't run Debug bconds, var_dt and absorber tests on MPI as it takes too long
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C kahan_sum test || cat kahan_sum/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C cone_bugs test || cat cone_bugs/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C shallow_water test || cat shallow_water/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then make -C concurrent_1d test || cat concurrent_1d/Testing/Temporary/LastTest.log /; fi
if [[ $MPI != 'none' ]]; then OMP_NUM_THREADS=2 make -C git_revision test || cat git_revision/Testing/Temporary/LastTest.log /; fi

cmake -DCMAKE_BUILD_TYPE=Release ../
# on MPI run absorber in release mode, bconds and var_dt take too long (why? they're fast on cuda-k-4)
if [[ $MPI != 'none' ]]; then VERBOSE=1 make -C absorber; fi
if [[ $MPI != 'none' ]]; then make -C absorber test || cat absorber/Testing/Temporary/LastTest.log /; fi

if [[ $MPI == 'none' ]]; then VERBOSE=1 $make_j; fi
# excluding test_issue because it (sometimes) fails in Release mode
# for unknown reasons and it only seems to happen on Travis ...
# "/" intentional! (just to make cat exit with an error code)
if [[ $MPI == 'none' ]]; then OMP_NUM_THREADS=4 ctest -E test_issue || cat Testing/Temporary/LastTest.log /;fi
cd ../../..
set +e
