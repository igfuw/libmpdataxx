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
# TODO: should't this path be added somewhere now?

# Boost libraries
find_package(Boost COMPONENTS thread date_time system iostreams timer filesystem REQUIRED QUIET)

# Blitz
find_package(Blitz++ REQUIRED) 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${BLITZ_INCLUDE_DIR}") # TODO: same for boost and others???

# HDF5 libraries
find_package(HDF5 COMPONENTS CXX HL)

# macro to be used in the subdirectories
function(libmpdataxx_add_test test)
  add_executable(${test} ${test}.cpp)
  set_target_properties(${test} PROPERTIES LINK_FLAGS_DEBUG "-lblitz") # TODO: use ${BLITZ_LIBRARIES}
  target_link_libraries(${test} ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
# TODO: target_compile_options() // added to CMake on Jun 3rd 2013
  add_test(${test} ${test})
endfunction()

# list of tests

#misc
add_subdirectory(kahan_sum)
# adv
add_subdirectory(paper_example_1)
add_subdirectory(paper_example_2)
add_subdirectory(gnuplot-iostream_1d)
add_subdirectory(gnuplot-iostream_2d)
add_subdirectory(var_sign_2d)
add_subdirectory(concurrent_1d)
add_subdirectory(err_isolines)
add_subdirectory(rotating_cone)
add_subdirectory(rotating_sphere)
add_subdirectory(mpdata_fct_1d)
add_subdirectory(mpdata_1d_opt_abs)
add_subdirectory(mpdata_1d_opt_tot)
add_subdirectory(mpdata_1d_opt_eps)
add_subdirectory(mpdata_1d_opt_npa)
add_subdirectory(mpdata_1d_opt_iga)
add_subdirectory(mpdata_1d_opt_nug)
# adv+rhs
add_subdirectory(harmosc)
# adv+rhs+vip
add_subdirectory(shallow_water)
# adv+rhs+vip+prs
add_subdirectory(bcbf)
