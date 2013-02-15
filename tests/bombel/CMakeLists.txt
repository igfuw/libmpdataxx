#TODO parallelisation does not work yet...
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_package(Boost COMPONENTS iostreams system REQUIRED)
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)

set(tmp test_bombel)
add_executable(${tmp} ${tmp}.cpp)
set_target_properties(${tmp} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")
target_link_libraries(${tmp} ${Boost_LIBRARIES})
add_test(${tmp} ${tmp})
