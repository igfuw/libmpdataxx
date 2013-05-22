cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(libmpdata++ CXX C)

# library
set(CMAKE_CXX_FLAGS "-std=c++11")
add_subdirectory(libmpdata++) 
#TODO: search for Blitz

# docs
find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_subdirectory(doc)
else()
  message(WARNING "Doxygent not found - 'make doc' will not be available.")
endif()

# tests
enable_testing()
add_subdirectory(tests)
