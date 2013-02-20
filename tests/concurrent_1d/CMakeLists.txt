set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#find_package(Threads REQUIRED) # does not work for clang :(
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # TODO (report to CMake?)

set(tmp test_concurrent_1d)
add_executable(${tmp} ${tmp}.cpp)
set_target_properties(${tmp} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")
target_link_libraries(${tmp} ${BOOST_THREAD_LIBS})
add_test(${tmp} ${tmp})
