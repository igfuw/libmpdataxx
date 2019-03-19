#!/usr/bin/env sh
set -e
cd tests/sandbox
mkdir build
cd build
cmake ..
VERBOSE=1 $make_j
if [[ $MPI != 'none' ]]; then  
  env OMP_NUM_THREADS=2
else
  env OMP_NUM_THREADS=4
fi

# running selected sandbox tests in Release mode
# "/" intentional! (just to make cat exit with an error code)
make -C mpi_adv test || cat mpi_adv/Testing/Temporary/LastTest.log /
ctest -R tgv_2d || cat Testing/Temporary/LastTest.log /
make -C convergence_2d_3d test || cat convergence_2d_3d/Testing/Temporary/LastTest.log /
travis_wait 30 make -C convergence_vip_1d test || cat convergence_vip_1d/Testing/Temporary/LastTest.log /
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then make -C convergence_spacetime test || cat convergence_spacetime/Testing/Temporary/LastTest.log /; fi
# with mpi it takes too long, bconds_div is ran in a separate test suite for mpi
if [[ $MPI != 'none' ]]; then make -C bconds_div test || cat bconds_div/Testing/Temporary/LastTest.log /; fi
cd ../../..
