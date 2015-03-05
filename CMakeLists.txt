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
  message(FATAL_ERROR "C++11 compatibility test failed - please update your compiler or point CMake to another one with -DCMAKE_CXX_COMPILER=...")
endif()

# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(COMMAND
   bash -c "git log -1 --format=\"format:#define GIT_REVISION \\\"%H\\\"%n\" HEAD > libmpdata++/git_revision.hpp"
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} 
  )
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
