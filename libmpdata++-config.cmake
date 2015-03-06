# TODO: optional components: hdf, gnuplot, openmp, threads, ...
#TODO: FindPackageHandleStandardArgs?


set(libmpdataxx_INCLUDE_DIRS "")
set(libmpdataxx_LIBRARIES "")
set(libmpdataxx_FOUND "")
set(libmpdataxx_CXX_FLAGS_DEBUG "")
set(libmpdataxx_CXX_FLAGS_RELEASE "")


############################################################################################
if(APPLE)
  # needed for the XCode clang to be identified as AppleClang and not Clang
  cmake_minimum_required(VERSION 3.0) 
  cmake_policy(SET CMP0025 NEW)
  cmake_policy(SET CMP0042 NEW)
else()
  # needed for the OpenMP test to work in C++-only project 
  # (see http://public.kitware.com/Bug/view.php?id=11910)
  cmake_minimum_required(VERSION 2.8.8) 
endif()


############################################################################################
# C++11
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-std=c++11")
check_cxx_source_compiles("
  #include <type_traits>
  template <bool a, class b> using ei=std::enable_if<a,b>; 
  struct a {a(int){}};struct b:a {using a::a;};  
  int main(){b i(1);}
" CXX11_SUPPORTED)
if (NOT CXX11_SUPPORTED)
  message(FATAL_ERROR "C++11 compatibility test failed - please update your compiler or point CMake to another one with -DCMAKE_CXX_COMPILER=...")
endif()


############################################################################################
# debug mode compiler flags
set(libmpdataxx_CXX_FLAGS_DEBUG "-std=c++11 -DBZ_DEBUG -g") #TODO: -Og if compiler supports it?


############################################################################################
# release mode compiler flags
if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR 
  CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
)
  set(libmpdataxx_CXX_FLAGS_RELEASE "-std=c++11 -DNDEBUG -Ofast -march=native")

  # preventing Kahan summation from being optimised out
  if (
    (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6) OR
    (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1) #TODO: never actually checked!
  )
    set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -fno-vectorize") 
  endif()
endif()


############################################################################################
# multi-threading
# find_package(ThreadsCXX) <- this requires C language to be enabled
# TODO: better solution!
set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -pthread")
set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -pthread")


############################################################################################
# Blitz++
find_path(BLITZ_INCLUDE_DIR blitz/array.h) #TODO: Fedora uses different path for config.h!
find_library(BLITZ_LIBRARY NAMES blitz)

include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_FLAGS "-I${BLITZ_INCLUDE_DIR}")
check_include_file_cxx("blitz/array.h" BLITZ_FOUND)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "Blitz++ includes not found - please install Blitz++ or point CMake to it.")
endif()

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <blitz/array.h>\nint main(){blitz::Array<float,1> a(2); a = blitz::safeToReturn(a+a);}" BLITZ_VER_AT_LEAST_0_10)
if (NOT BLITZ_VER_AT_LEAST_0_10)
  message(FATAL_ERROR "Blitz++ ver >= 0.10 requirement not met - please update Blitz++.")
endif()

set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${BLITZ_INCLUDE_DIR}")
set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${BLITZ_LIBRARY}")

############################################################################################
# OpenMP
find_package(OpenMP QUIET)
set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")


############################################################################################
# Boost libraries
find_package(Boost COMPONENTS thread date_time system iostreams timer filesystem QUIET)
if (NOT Boost_FOUND)
  message(FATAL_ERROR "
  Boost.{Thread,date_time,system,iostreams,timer,filesystem} not found.
  Please install it (e.g. sudo apt-get install libboost-all-dev).
")
endif()
set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${Boost_LIBRARIES}")
#TODO: include_dirs


############################################################################################
# HDF5 libraries
find_package(HDF5 COMPONENTS CXX HL) # REQUIRED?
set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${HDF5_LIBRARIES}")
#TODO: include_dirs


############################################################################################
# gnuplot-iostream
#TODO: check if gnuplot installed
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)
if (NOT GNUPLOT-IOSTREAM_INCLUDE_DIR)
  message(FATAL_ERROR "
  gnuplot-iostream not found. Please install it, e.g.:
    Debian/Ubuntu: sudo apt-get install libgnuplot-iostream-dev
    manual: wget -O /usr/local/include/gnuplot-iostream.h http://gitorious.org/gnuplot-iostream/gnuplot-iostream/raw/gnuplot-iostream.h 
")
endif()
set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${GNUPLOT-IOSTREAM_INCLUDE_DIR}")


############################################################################################
set(libmpdataxx_FOUND TRUE)
