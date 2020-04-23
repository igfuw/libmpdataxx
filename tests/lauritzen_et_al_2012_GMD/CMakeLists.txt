if(APPLE)
  # needed for the XCode clang to be identified as AppleClang and not Clang
  cmake_minimum_required(VERSION 3.0) 
else()
  # needed for the OpenMP test to work in C++-only project 
  # (see http://public.kitware.com/Bug/view.php?id=11910)
  cmake_minimum_required(VERSION 2.8.8) 
endif()

project(libmpdata++-tests-lauritzen_et_al_2012 CXX)

include(${CMAKE_SOURCE_DIR}/../../libmpdata++-config.cmake)
if(NOT libmpdataxx_FOUND) 
  message(FATAL_ERROR "local libmpdata++-config.cmake not found!")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libmpdataxx_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG ${libmpdataxx_CXX_FLAGS_DEBUG})
endif()

# to make <libmpdata++/...> work
set(CMAKE_CXX_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CXX_FLAGS}")

# macro to be used in the subdirectories
function(libmpdataxx_add_test test)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${libmpdataxx_LIBRARIES})
  target_include_directories(${test} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
  if(USE_MPI)
    add_test(NAME ${test} COMMAND ${libmpdataxx_MPIRUN} -np 3 ${CMAKE_CURRENT_BINARY_DIR}/${test}) # TODO: why 3? :)
  else()
    add_test(${test} ${test})
  endif()
endfunction()

enable_testing()

libmpdataxx_add_test(gmd_2012)

add_test(gmd_2012_stats_and_plots bash -c "
  python  ${CMAKE_CURRENT_SOURCE_DIR}/stats_and_plots.py nug_i2_120 nug_i2_240 nug_iga_fct_i2_120 nug_iga_fct_i2_240;
") 

add_test(gmd_2012_stats_diff bash -c "
  for i in nug_i2_120 nug_i2_240 nug_iga_fct_i2_120 nug_iga_fct_i2_240; do 
    echo  ${CMAKE_CURRENT_SOURCE_DIR}/stats_$i.txt.gz;
    echo  ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt;
    zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats_$i.txt.gz ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt || exit 1; 
  done
") 
