cmake_minimum_required(VERSION 2.8.12)
project(libmpdata++ CXX)

# using include() istead of find_package(libmpdata++) to use local CMake code
# and not the system-installed one
include(${CMAKE_SOURCE_DIR}/../libmpdata++-config.cmake)
if (NOT libmpdataxx_FOUND) 
  message(FATAL_ERROR "local libmpdata++-config.cmake not found!")
endif()

# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/../.git")
  execute_process(COMMAND
   bash -c "git log -1 --format=\"format:#define GIT_REVISION \\\"%H\\\"%n\" HEAD > git_revision.hpp"
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} 
  )
endif()

install(
  DIRECTORY
    bcond concurr formulae output solvers
  DESTINATION 
    include/libmpdata++
)
install(
  FILES
    blitz.hpp git_revision.hpp kahan_reduction.hpp opts.hpp
  DESTINATION 
    include/libmpdata++
)
install(
  FILES
    ../libmpdata++-config.cmake
  DESTINATION
    share/libmpdata++
)
