cmake_minimum_required(VERSION 3.1)

project(libmpdata++-tests-mp3_paper_2018_JCP CXX)

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

if (NOT FULL_SIM)
  set(FULL_SIM "false")
endif()

# macro to be used in the subdirectories
function(libmpdataxx_add_test test)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${libmpdataxx_LIBRARIES})
  target_include_directories(${test} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
  target_compile_definitions(${test} PRIVATE FULL_SIM=${FULL_SIM})
  add_test(${test} ${test})
endfunction()

enable_testing()

add_subdirectory(manufactured_3d)
add_subdirectory(reversing_deform)
add_subdirectory(moving_vort)
