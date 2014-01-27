cmake_minimum_required(VERSION 2.8.8) # for OpenMP test to work in C++-only project (see http://public.kitware.com/Bug/view.php?id=11910)
project(libmpdata++ CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# C++11
set(CMAKE_CXX_FLAGS "-std=c++11")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("struct a {a(int){}};struct b:a {using a::a;};int main(){b i(1);}" CXX11_SUPPORTED)
if (NOT CXX11_SUPPORTED)
  message(FATAL_ERROR "C++11 compatibility test failed - please update your compiler.")
endif()

# Blitz
find_package(Blitz++ REQUIRED) #TODO: should it be here or in tests?

# the library
add_subdirectory(libmpdata++) 

# the docs
find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_subdirectory(doc)
else()
  message(WARNING "Doxygent not found - 'make doc' will not be available.")
endif()

# the tests
enable_testing()
add_subdirectory(tests)
