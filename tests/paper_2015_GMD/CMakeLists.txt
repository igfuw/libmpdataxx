cmake_minimum_required(VERSION 2.8.12)
project(libmpdata++-tests-paper_2015_GMD CXX)

# using include() istead of find_package(libmpdata++) to use local CMake code
# and not the system-installed one
include(${CMAKE_SOURCE_DIR}/../../libmpdata++-config.cmake)
if(NOT libmpdataxx_FOUND) 
  message(FATAL_ERROR "local libmpdata++-config.cmake not found!")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS ${libmpdataxx_CXX_FLAGS_RELEASE})
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG ${libmpdataxx_CXX_FLAGS_DEBUG})
endif()

# to make <libmpdata++/...> work
set(CMAKE_CXX_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CXX_FLAGS}")

# TODO: check for Paraview (+ version requirement)
# TODO: check for matplotlib, ...

# gzip
execute_process(
  COMMAND bash -c "
    TMP=`mktemp cmake-zdiff-XXXXX` &&
    gzip < $TMP > $TMP.gz          &&
    zdiff $TMP.gz $TMP             &&
    rm $TMP $TMP.gz
  " 
  OUTPUT_QUIET ERROR_QUIET
  RESULT_VARIABLE STATUS
)
if(NOT STATUS EQUAL 0)
  SET(msg "zdiff does not work properly!")
  if(APPLE)
    message(FATAL_ERROR "${msg}
      On OSX, if using Homebrew, the foloowing should help:
        brew tap homebrew/dupes
        brew install gzip
        export PATH=/usr/local/bin:$PATH
    ")
  else()
    message(FATAL_ERROR ${msg})
  endif()
endif()

# macro to be used in the subdirectories
function(libmpdataxx_add_test test)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${libmpdataxx_LIBRARIES})
  target_include_directories(${test} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
  add_test(${test} ${test})
endfunction()

function(libmpdataxx_add_test_gi test)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${libmpdataxx_LIBRARIES})
  target_include_directories(${test} PUBLIC ${libmpdataxx_INCLUDE_DIRS})

  # this will not work with multiple threads (make -j) as it depends on the order of execution + $PPID might be different?
  add_test(NAME "init_${test}" COMMAND "bash" "-c" "> log-$PPID")
  add_test(NAME "calc_${test}" COMMAND "bash" "-c" "GNUPLOT_IOSTREAM_CMD=\"cat >> log-`echo $PPID`\" ${CMAKE_CURRENT_BINARY_DIR}/${test}")
  add_test(NAME "plot_${test}" COMMAND "bash" "-c" "gnuplot log-$PPID")
  add_test(NAME "diff_${test}" COMMAND "bash" "-c" "zdiff log-$PPID ${CMAKE_CURRENT_SOURCE_DIR}/refdata/log-orig.gz && rm log-$PPID")
endfunction()

enable_testing()

# adv
add_subdirectory(0_basic_example)
add_subdirectory(1_advscheme_opts)
add_subdirectory(2_convergence_1d)
add_subdirectory(3_rotating_cone_2d)
add_subdirectory(4_revolving_sphere_3d)
add_subdirectory(5_over_the_pole_2d)
# adv+rhs
add_subdirectory(6_coupled_harmosc)
# adv+rhs+vip
add_subdirectory(7_shallow_water)
# adv+rhs+vip+prs
add_subdirectory(8_boussinesq_2d)
