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
if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR 
  CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"
)
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -DNDEBUG -Ofast -march=native")

  # preventing Kahan summation from being optimised out
  if (
    (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5) OR
    (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1) #TODO: never actually checked!
  )
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-vectorize") 
  endif()
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
find_package(HDF5 COMPONENTS CXX HL REQUIRED)

# macro to be used in the subdirectories
function(libmpdataxx_add_test test)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
  target_link_libraries(${test} debug "-lblitz") # TODO: use ${BLITZ_LIBRARIES}
  target_include_directories(${test} PUBLIC ${HDF5_INCLUDE_DIRS})
# TODO: target_compile_options() // added to CMake on Jun 3rd 2013
  add_test(${test} ${test})
endfunction()

function(libmpdataxx_add_test_gi test)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${Boost_LIBRARIES})
  target_link_libraries(${test} debug "-lblitz") # TODO: use ${BLITZ_LIBRARIES}

  # this will not work with multiple threads (make -j) as it depends on the order of execution + $PPID might be different?
  add_test(NAME "init_${test}" COMMAND "bash" "-c" "> log-$PPID")
  add_test(NAME "calc_${test}" COMMAND "bash" "-c" "GNUPLOT_IOSTREAM_CMD=\"cat >> log-`echo $PPID`\" ${CMAKE_CURRENT_BINARY_DIR}/${test}")
  add_test(NAME "plot_${test}" COMMAND "bash" "-c" "gnuplot log-$PPID")
  add_test(NAME "diff_${test}" COMMAND "bash" "-c" "zcmp log-$PPID ${CMAKE_CURRENT_SOURCE_DIR}/log-orig.gz && rm log-$PPID")
endfunction()

# tests
add_subdirectory(tutorial)
add_subdirectory(sandbox)
add_subdirectory(unit)
