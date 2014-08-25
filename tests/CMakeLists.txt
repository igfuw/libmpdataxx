# to make <libmpdata++/...> work
set(CMAKE_CXX_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CXX_FLAGS}") 

# OpenMP
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#find_package(Threads REQUIRED) # does not work for clang :(
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # TODO (report to CMake?)

# Debug mode
set(CMAKE_CXX_FLAGS_DEBUG "-DBZ_DEBUG -g")

# Release mode
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -DNDEBUG -Ofast -march=native")
endif()

# Profile mode
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -pg")

# gnuplot-iostream
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)
if (NOT GNUPLOT-IOSTREAM_INCLUDE_DIR)
  message(FATAL_ERROR "
  gnuplot-iostream not found.
  Please install it (e.g. sudo apt-get install libgnuplot-iostream-dev).
")
endif()
# TODO: should't this path be added somewhere now?

# Boost libraries
find_package(Boost COMPONENTS thread date_time system iostreams timer filesystem QUIET)
if (NOT Boost_FOUND)
  message(FATAL_ERROR "
  Boost.{Thread,date_time,system,iostreams,timer,filesystem} not found.
  Please install it (e.g. sudo apt-get install libboost-all-dev).
")
endif()

# Blitz
find_package(Blitz++ REQUIRED) 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${BLITZ_INCLUDE_DIR}") # TODO: same for boost and others???

# HDF5 libraries
find_package(HDF5 COMPONENTS CXX HL)

# macro to be used in the subdirectories
function(libmpdataxx_add_test test)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
  target_link_libraries(${test} debug "-lblitz") # TODO: use ${BLITZ_LIBRARIES}
# TODO: target_compile_options() // added to CMake on Jun 3rd 2013
  add_test(${test} ${test})
endfunction()

# tests
add_subdirectory(tutorial)
add_subdirectory(sandbox)
add_subdirectory(unit)
