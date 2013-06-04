# to make <libmpdata++/...> work
set(CMAKE_CXX_FLAGS "-I../.. ${CMAKE_CXX_FLAGS}") 

# OpenMP
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#find_package(Threads REQUIRED) # does not work for clang :(
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # TODO (report to CMake?)

# Blitz Debug mode
set(CMAKE_CXX_FLAGS_DEBUG "-DBZ_DEBUG -g")

# -DNDEBUG?
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -Ofast -march=native")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -O3 -mllvm -vectorize -ffast-math -march=native") 
endif()

# gnuplot-iostream
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)
# TODO: should't this path be added somewhere now?

# Boost libraries
find_package(Boost COMPONENTS thread date_time system iostreams timer REQUIRED QUIET)

# macro to be used in the subdirectories
function(libmpdataxx_add_test test)
  add_executable(${test} ${test}.cpp)
  set_target_properties(${test} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")
  target_link_libraries(${test} ${Boost_LIBRARIES})
  add_test(${test} ${test})
endfunction()

# list of tests
add_subdirectory(kahan_sum)
add_subdirectory(gnuplot-iostream_1d)
add_subdirectory(gnuplot-iostream_2d)
add_subdirectory(var_sign_2d)
add_subdirectory(harmosc)
add_subdirectory(todo_3d)
add_subdirectory(bombel)
add_subdirectory(concurrent_1d)
add_subdirectory(shallow_water)
add_subdirectory(rotating_cone)
#add_subdirectory(bcbf)
add_subdirectory(mpdata_fct_1d)
add_subdirectory(mpdata_1d_opt_sss)
add_subdirectory(mpdata_1d_opt_toa)
