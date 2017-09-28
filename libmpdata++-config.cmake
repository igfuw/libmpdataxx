if(APPLE)
  # needed for the XCode clang to be identified as AppleClang and not Clang
  cmake_minimum_required(VERSION 3.0) 
else()
  # needed for the OpenMP test to work in C++-only project 
  # (see http://public.kitware.com/Bug/view.php?id=11910)
  cmake_minimum_required(VERSION 2.8.8) 
endif()

# the policies we care about:
# - CMP0025 - make CMake distinguis between Apple and LLVM clang
# - CMP0042 - make CMake use RPATHs on OSX
if(CMAKE_VERSION VERSION_GREATER 2.9)
  cmake_policy(VERSION 3.0)
endif()

############################################################################################
# the following variables will be set:
set(libmpdataxx_FOUND False)
set(libmpdataxx_INCLUDE_DIRS "")
set(libmpdataxx_LIBRARIES "")
set(libmpdataxx_CXX_FLAGS_DEBUG "")
set(libmpdataxx_CXX_FLAGS_RELEASE "")

############################################################################################
# libmpdata++ headers for non-default install location (i.e. for make DESTDIR=<dir> install)
set(libmpdataxx_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../include/")


############################################################################################
# debug mode compiler flags
set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -std=c++14 -DBZ_DEBUG -g") #TODO: -Og if compiler supports it?


############################################################################################
# release mode compiler flags
if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR 
  CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
)
  set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -std=c++14 -DNDEBUG -Ofast -march=native")

  # preventing Kahan summation from being optimised out
  if (
    (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6) OR
    (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1) #TODO: never actually checked!
  )
    set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -fno-vectorize") 
  endif()
endif()


############################################################################################
# C++14
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-std=c++14")
check_cxx_source_compiles("
  #include <type_traits>
  auto f() { return 1;}
  template <bool a, class b> using ei=std::enable_if<a,b>; 
  struct a {a(int){}};struct b:a {using a::a;};  
  int main(){b i(1);}
" CXX14_SUPPORTED)
if (NOT CXX14_SUPPORTED)
  message(FATAL_ERROR "C++14 compatibility test failed - please update your compiler or point CMake to another one with -DCMAKE_CXX_COMPILER=...")
endif()
unset(CMAKE_REQUIRED_FLAGS)


############################################################################################
# Blitz++
find_package(PkgConfig)
pkg_check_modules(BLITZ QUIET blitz>=0.10)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "Blitz++ library not found or Blitz++ version requirement not met (>=0.10)

* To insall Blitz++, please try:
*   Debian/Ubuntu: sudo apt-get install libblitz0-dev
*   Fedora: sudo yum install blitz-devel
*   Homebrew: brew install blitz
  ")
else()
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${BLITZ_INCLUDE_DIRS}")
  set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${BLITZ_LIBRARIES}")
endif()

############################################################################################
# OpenMP
find_package(OpenMP QUIET)
if(OPENMP_FOUND)
  set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
  set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
else()
  message(STATUS "OpenMP not supported by the compiler.

* Programs using libmpdata++'s OpenMP concurrency will compile but run in serial mode
  ")
endif()


############################################################################################
# multi-threading
# find_package(ThreadsCXX) <- this requires C language to be enabled
# TODO: better solution! 
# TODO: not needed for serial-only programs!
# TODO: -fopenmp implies -pthread on gcc
set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -pthread")
set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -pthread")


############################################################################################
# Boost libraries v>=1.55.0, because boost/predef was added then
set(Boost_DETAILED_FAILURE_MSG ON)
find_package(Boost 1.55.0 COMPONENTS thread date_time system iostreams timer filesystem)
if(Boost_FOUND)
  set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${Boost_LIBRARIES}")
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${Boost_INCLUDE_DIRS}")
  if(
    (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") 
    AND (Boost_MINOR_VERSION EQUAL 55) 
  )
    # add a definition -DBOOST_HAS_INT128=1 to clang calls on linux to avoid errors with boost.atomic (https://svn.boost.org/trac/boost/ticket/9610)
    set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -DBOOST_HAS_INT128=1")
    set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -DBOOST_HAS_INT128=1")
  endif()
else()
  #TODO: check separately for optional and mandatory components
  message(FATAL_ERROR "Boost (or some of its components) not found.

* Programs based on libmpdata++ will not compile. 
* To insall Boost, please try:
*   Debian/Ubuntu: sudo apt-get install libboost-all-dev
*   Fedora: sudo yum install boost-devel
  ")
endif()


############################################################################################
# HDF5 libraries
find_package(HDF5 COMPONENTS CXX HL QUIET)
if(HDF5_FOUND)
  set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${HDF5_LIBRARIES}")
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${HDF5_INCLUDE_DIRS}")
else()
  message(STATUS "HDF5 not found. 

* Programs using libmpdata++'s HDF5 output will not compile.
* To install HDF5, please try:
*   Debian/Ubuntu: sudo apt-get install libhdf5-serial-dev hdf5-tools
*   Fedora: sudo yum install hdf5-devel
*   Homebrew: brew install hdf5 --with-cxx
  ")
endif()


############################################################################################
# gnuplot-iostream
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)
if(GNUPLOT-IOSTREAM_INCLUDE_DIR)
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${GNUPLOT-IOSTREAM_INCLUDE_DIR}")
else()
  message(STATUS "gnuplot-iostream not found.

* Programs using libmpdata++'s gnuplot-iostream output will not compile.
* To install gnuplot-iostream, please try:
*   Debian/Ubuntu: sudo apt-get install libgnuplot-iostream-dev
*   manual: wget -O /usr/local/include/gnuplot-iostream.h https://raw.githubusercontent.com/dstahlke/gnuplot-iostream/master/gnuplot-iostream.h
  ")
endif()

find_program(GNUPLOT_FOUND NAMES gnuplot)
if(GNUPLOT_FOUND)
else()
  message(STATUS "gnuplot not found.

* Programs using libmpdata++'s gnuplot-iostream output will not run.
* To install gnuplot, please try:
*   Debian/Ubuntu: sudo apt-get install gnuplot
*   Fedora: sudo yum install gnuplot
*   Homebrew: brew install gnuplot
  ")
endif()


############################################################################################
list(REMOVE_DUPLICATES libmpdataxx_INCLUDE_DIRS)
list(REMOVE_ITEM libmpdataxx_INCLUDE_DIRS "")
set(libmpdataxx_FOUND TRUE)
