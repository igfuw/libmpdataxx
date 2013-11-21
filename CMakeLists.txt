cmake_minimum_required(VERSION 2.8.0) # TODO: version requirement took from another context!
project(libmpdata++ CXX)

# C++11
set(CMAKE_CXX_FLAGS "-std=c++11")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("struct a {a(int){}};struct b:a {using a::a;};int main(){b i(1);}" CXX11_SUPPORTED)
if (NOT CXX11_SUPPORTED)
  message(FATAL_ERROR "C++11 compatibility test failed - please update your compiler.")
endif()

# Blitz
include(CheckIncludeFileCXX)
check_include_file_cxx("blitz/array.h" BLITZ_FOUND)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "Blitz++ includes not found - please install Blitz++ or point CMake to it.")
endif()
check_cxx_source_compiles("#include <blitz/array.h>\nint main(){blitz::Array<float,1> a(2); a = blitz::safeToReturn(a+a);}" BLITZ_VER_AT_LEAST_0_10)
if (NOT BLITZ_VER_AT_LEAST_0_10)
  message(FATAL_ERROR "Blitz++ ver >= 0.10 requirement not met - please update Blitz++.")
endif()

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
