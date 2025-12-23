if(APPLE)
  # needed for the XCode clang to be identified as AppleClang and not Clang
  cmake_minimum_required(VERSION 3.0) 
else()
  # needed for the OpenMP test to work in C++-only project 
  # (see http://public.kitware.com/Bug/view.php?id=11910)
  cmake_minimum_required(VERSION 2.8.8) 
endif()

project(libmpdata++-tests-unit CXX)

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
#  message("libmpdataxx_LIBRARIES: ${libmpdataxx_LIBRARIES}")
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

# add_subdirectory(kahan_sum) # disabled, because since uwlcm_ubuntu_24_04_cuda_12_9_0.sif, blitz::sum gives the same results as kahan sum; failed e.g. in uwlcm_ubuntu_20_04_cuda_11_4.sif
add_subdirectory(cone_bugs)
add_subdirectory(shallow_water)
add_subdirectory(concurrent_1d)
if(!USE_MPI)
  add_subdirectory(hint_scale) # initialization from pre-defined arrays, not using index placeholders
  add_subdirectory(hdf5_catch) # parallel_hdf5 exceptions - couldn't find any documentation
endif()
add_subdirectory(git_revision)
add_subdirectory(absorber)
add_subdirectory(bconds)
add_subdirectory(var_dt)
add_subdirectory(delayed_advection)
