cmake_minimum_required(VERSION 2.8.8) # for OpenMP test to work in C++-only project (see http://public.kitware.com/Bug/view.php?id=11910)
project(libmpdata++ CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
  #include <type_traits>
  template <bool a, class b> using ei=std::enable_if<a,b>; 
  struct a {a(int){}};struct b:a {using a::a;};  
  int main(){b i(1);}
" CXX11_SUPPORTED)
if (NOT CXX11_SUPPORTED)
  message(FATAL_ERROR "C++11 compatibility test failed - please update your compiler.")
endif()

# the library
add_subdirectory(libmpdata++) 

# the tests
enable_testing()
add_subdirectory(tests)

# make-install rule #TODO: exclude CmakeLists.txt
install(
  DIRECTORY
    libmpdata++
  DESTINATION 
    include
)

# make-install rule for FindLibmpdataxx.cmake
install(
  FILES
    libmpdata++-config.cmake
  DESTINATION
    share/libmpdata++
)
