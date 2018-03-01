
libmpdata++ - a library of parallel MPDATA-based solvers for systems of generalised transport equations 
=======================================================================

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c1be1645743a463b8ead70f111d2ba15)](https://app.codacy.com/app/igfuw/libmpdataxx?utm_source=github.com&utm_medium=referral&utm_content=igfuw/libmpdataxx&utm_campaign=badger)

To get more information on libmpdata++, please check: 
  - http://libmpdataxx.igf.fuw.edu.pl/
  - http://arxiv.org/abs/1407.1309
  - http://www.geosci-model-dev.net/8/1005/2015/

In short, libmpdata++ is a header-only C++ library. 
Compilation of programs that use libmpdata++ requires:
- a C++11 compliant compiler (optionally with OpenMP support)
- Blitz++ and Boost C++ libraries
- HDF5 and gnuplot-iostream libraries 
  (optional, depending on the type of output mechanism chosen)

During development of libmpdata++, we are continuously testing
the code on Linux using GCC and LLVM/Clang as well as on OSX
using Apple/Clang - these are considered the supported platforms.

Compilation and execution of the examples shipped with libmpdata++ 
is easiest done using CMake, and the following instructions assume
you're using CMake. Some hints on CMake usage are included at the
end of this file.

The .travis.yml file shipped with the library contains a complete
set of commands needed to build and execute all tests programs
shipped with libmpdata++ on fresh Ubuntu and OSX installations -
it may contain useful information on obtaining the dependencies.

1. To verify if all dependencies are met, please start with:
  $ cd libmpdata++
  $ mkdir build
  $ cd build
  $ cmake ..
  $ cd ../..
  
The next two steps are optional test. Running the tests is highly
recommended to verify if the library works correctly in your 
environment. Nevertheless, in principle you can skip to step four
and install the library right away.
  
2. To perform unit tests, please try:
  $ cd tests/unit
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make
  $ make test
  $ cd ../../..

The unit tests should complete in a dozen of seconds.

3. To reproduce all results from the GMD paper, please try:
  $ cd tests/paper_2015_GMD
  $ mkdir build 
  $ cd build
  $ cmake ..
  $ make
  $ make test     
  $ cd ../../..

This takes ca. 15 minutes on a quad-core laptop. The "make test"
command performs simulations, checks the output against reference 
data (tests/paper_2015_GMD/*/refdata/*) and plots all figures 
included in the paper. The subfolders of paper_2015_GMD correspond 
to consecutive chapters in the GMD paper. Some of the scripts run
by "make test" require additional packages including Python, Python
libraries (NumPy, SciPy, matplotlib) and Paraview.

4. To install the library system-wide, please try:
  $ cd libmpdata++/build
  $ sudo make install

This will copy the libmpdata++ headers into the system include path
(e.g. /usr/include/libmpdata++) and copy the libmpdata++-config.cmake 
file into the system share directory (e.g. /usr/share/libmpdata++) 
what will allow CMake users to do find_package(libmpdata++).

Some CMake hints:
- to point CMake to a non-default C++ compiler (e.g. clang++):
  $ cmake .. -DCMAKE_CXX_COMPILER=clang++ 

- to alter the installation prefix (e.g. /usr/ instead of /usr/local):
  $ cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/usr

- to switch between debug and release (default) compilation modes 
  (has to be done after compiler choice):
  $ cmake .. -DCMAKE_BUILD_TYPE=Debug
  $ cmake .. -DCMAKE_BUILD_TYPE=Release
  
- two alternative ways of cleaning leftovers from a previous build 
  (including CMake cache files):
  $ rm -rf build/CMakeCache.txt build/CMakeFiles
  $ rm -rf build; mkdir build

- the output of commands executed by "make test" can be viewed with:
  $ less Testing/Temporary/LastTest.log
