set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_package(Boost COMPONENTS thread date_time system REQUIRED)
find_package(Threads REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # TODO (report to CMake?)

set(tmp test_concurrent_1d)
add_executable(${tmp} ${tmp}.cpp)
set_target_properties(${tmp} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")
#target_link_libraries(${tmp} ${CMAKE_THREAD_LIBS_INIT}) # see above
target_link_libraries(${tmp} ${Boost_LIBRARIES})
add_test(${tmp} ${tmp})
