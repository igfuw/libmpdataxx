set(CMAKE_CXX_FLAGS "-I../.. ${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_DEBUG "-DBZ_DEBUG -Wfatal-errors")

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Ofast -march=native")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")
endif()

find_package(OpenMP)

add_subdirectory(gnuplot-iostream_1d)
add_subdirectory(gnuplot-iostream_2d)
add_subdirectory(var_sign_2d)
add_subdirectory(harmosc)
add_subdirectory(todo_3d)
add_subdirectory(bombel)
add_subdirectory(concurrent_1d)
