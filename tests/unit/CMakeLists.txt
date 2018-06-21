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
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${libmpdataxx_LIBRARIES})
  target_include_directories(${test} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
  add_test(${test} ${test})
endfunction()

enable_testing()

add_subdirectory(kahan_sum)
add_subdirectory(cone_bugs)
add_subdirectory(shallow_water)
add_subdirectory(concurrent_1d)
add_subdirectory(hint_scale)
add_subdirectory(hdf5_catch)
add_subdirectory(git_revision)
add_subdirectory(absorber)
add_subdirectory(bconds)
add_subdirectory(var_dt)
add_subdirectory(delayed_advection)
