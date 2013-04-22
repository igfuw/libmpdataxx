set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

find_package(Boost COMPONENTS iostreams system timer REQUIRED)
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I../../../libcloudphxx") # TODO: handle it correctly!!!

set(tmp test_lcph_kessler)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} ${Boost_LIBRARIES})
set_target_properties(${tmp} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")
target_link_libraries(${tmp} ${BOOST_THREAD_LIBS})
add_test(${tmp} ${tmp})
