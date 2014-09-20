include(CheckCXXSourceCompiles)

# C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
check_cxx_source_compiles("
  #include <type_traits>
  template <bool a, class b> using ei=std::enable_if<a,b>; 
  struct a {a(int){}};struct b:a {using a::a;};  
  int main(){b i(1);}
" CXX11_SUPPORTED)
if (NOT CXX11_SUPPORTED)
  message(FATAL_ERROR "C++11 compatibility test failed - please update your compiler.")
endif()

# pthreads (TODO: more portable way of setting it)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

# Debug mode
set(CMAKE_CXX_FLAGS_DEBUG "-DBZ_DEBUG -g")

# Release mode
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -DNDEBUG -Ofast -march=native")
endif()

# Boost libraries
find_package(Boost COMPONENTS thread date_time system iostreams timer filesystem QUIET)
if (NOT Boost_FOUND)
  message(FATAL_ERROR "
  Boost.{Thread,date_time,system,iostreams,timer,filesystem} not found.
  Please install it (e.g. sudo apt-get install libboost-all-dev).
")
endif()

# HDF5 libraries
find_package(HDF5 COMPONENTS CXX HL)

set(libmpdataxx_LIBRARIES "${Boost_LIBRARIES};${HDF5_LIBRARIES}")
message("${libmpdataxx_LIBRARIES}")
set(libmpdataxx_FOUND TRUE)
