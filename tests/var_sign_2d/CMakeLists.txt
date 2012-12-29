find_package(Boost COMPONENTS iostreams system REQUIRED)
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)

set(tmp test_var_sign_2d)
add_executable(${tmp} ${tmp}.cpp)
target_link_libraries(${tmp} ${Boost_LIBRARIES})
add_test(${tmp} ${tmp})
