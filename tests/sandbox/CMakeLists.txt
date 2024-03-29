if(APPLE)
  # needed for the XCode clang to be identified as AppleClang and not Clang
  cmake_minimum_required(VERSION 3.0) 
else()
  # needed for the OpenMP test to work in C++-only project 
  # (see http://public.kitware.com/Bug/view.php?id=11910)
  cmake_minimum_required(VERSION 2.8.8) 
endif()

project(libmpdata++-tests-sandbox CXX)

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

# macros to be used in the subdirectories
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

add_subdirectory(mpi_adv)
add_subdirectory(straka)
add_subdirectory(tgv)
add_subdirectory(convergence_2d_3d)
add_subdirectory(pbl)
add_subdirectory(convergence_spacetime)
add_subdirectory(bconds_div)
add_subdirectory(shear_layer)
add_subdirectory(convergence_vip_1d)
add_subdirectory(convergence_adv_diffusion)
