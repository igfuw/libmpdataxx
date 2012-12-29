cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(advoocat CXX)

# library
set(CMAKE_CXX_FLAGS "-std=c++0x -Wfatal-errors")
add_subdirectory(advoocat) 

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
