libmpdata++ - a library of parallel MPDATA-based solvers for systems of generalised transport equations
=======================================================================

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c1be1645743a463b8ead70f111d2ba15)](https://app.codacy.com/app/igfuw/libmpdataxx?utm_source=github.com&utm_medium=referral&utm_content=igfuw/libmpdataxx&utm_campaign=badger)

MPDATA stands for Multidimensional Positive Definite Advection Transport Algorithm.
It applies to a variety of problems involving conservation laws in computational fluid dynamics.
More generally, it can be used for numerically tackling linear hyperbolic first-order
partial differential equation systems. The algorithm was introduced in the context of simulations of atmospheric flows
([Smolarkiewicz, 1984](http://doi.org/10.1016/0021-9991(84)90121-9)).

MPDATA is explicit, forward-in-time, sign-preserving, conservative and non-linearly stable.
Through iterative corrections, it achieves high-order accuracy in time and space even for complex flows and problem geometries.
The algorithm has been continuously developed over a third of a century
resulting in a family of robust numerical schemes with documented applications in computational bio-, geo-, and astro-physics as well as engineering
(see [Smolarkiewicz, 2006](http://doi.org/10.1002/fld.1071), [Smolarkiewicz et al., 2016](http://doi.org/10.1016/j.jcp.2016.06.048),
for an overview and recent references).

In short, libmpdata++ is a new header-only C++ implementation of MPDATA.
The features of the library along with discussion of several example
use cases were published in a user-guide-style journal paper:
[Jaruga et al., 2015](http://www.geosci-model-dev.net/8/1005/2015/).

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
We also provide a Docker image that contains all requirements
and the libmpdata++ library.

## Installing the library and running tests locally

### 1. To verify if all dependencies are met, please start with:
```
  $ cd libmpdata++
  $ mkdir build
  $ cd build
  $ cmake ..
  $ cd ../..
```
The next two steps are optional test. Running the tests is highly
recommended to verify if the library works correctly in your
environment. Nevertheless, in principle you can skip to step four
and install the library right away.

### 2. To perform unit tests:
```
  $ cd tests/unit
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make
  $ make test
  $ cd ../../..
```
The unit tests should complete in a dozen of seconds.

### 3. To reproduce all results from the GMD paper:
```
  $ cd tests/paper_2015_GMD
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make
  $ make test
  $ cd ../../..
```
This takes ca. 15 minutes on a quad-core laptop. The "make test"
command performs simulations, checks the output against reference
data (tests/paper_2015_GMD/*/refdata/*) and plots all figures
included in the paper. The subfolders of paper_2015_GMD correspond
to consecutive chapters in the GMD paper. Some of the scripts run
by "make test" require additional packages including Python, Python
libraries (NumPy, SciPy, matplotlib) and Paraview.

### 4. To install the library system-wide, please try:
```
  $ cd libmpdata++/build
  $ sudo make install
```
This will copy the libmpdata++ headers into the system include path
(e.g. /usr/include/libmpdata++) and copy the libmpdata++-config.cmake
file into the system share directory (e.g. /usr/share/libmpdata++)
what will allow CMake users to do find_package(libmpdata++).

### Some CMake hints:
- to point CMake to a non-default C++ compiler (e.g. clang++):
```
  $ cmake .. -DCMAKE_CXX_COMPILER=clang++
```
- to alter the installation prefix (e.g. /usr/ instead of /usr/local):
```
  $ cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/usr
```
- to switch between debug and release (default) compilation modes
  (has to be done after compiler choice):
```
  $ cmake .. -DCMAKE_BUILD_TYPE=Debug
  $ cmake .. -DCMAKE_BUILD_TYPE=Release
```
- two alternative ways of cleaning leftovers from a previous build
  (including CMake cache files):
```
  $ rm -rf build/CMakeCache.txt build/CMakeFiles
  $ rm -rf build; mkdir build
```
- the output of commands executed by "make test" can be viewed with:
```
  $ less Testing/Temporary/LastTest.log
```

## Using and testing the library with a [Docker](https://docs.docker.com/) image

- In order to use the image you have to [install Docker](https://docs.docker.com/install/).
Once you have Docker, the image can be downloaded:

```
docker pull igfuw/libmpdataxx:latest
```
You can also download an image with a specific version of the library,
a full list of available tags can be found [here](https://hub.docker.com/r/igfuw/libmpdataxx/tags/).

- To run the Docker in interactive mode:
```
docker run -it --rm igfuw/libmpdataxx:latest
```
This will open an interactive `bash` in the container
and your working directory will be `/usr/local/src/libmpdataxx`.
You can repeat the [step 2](#2-to-perform-unit-tests)
and the [step 3](#3-to-reproduce-all-results-from-the-gmd-paper) from the previous part to run the tests.
