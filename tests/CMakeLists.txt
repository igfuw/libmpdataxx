set(CMAKE_CXX_FLAGS "-I../.. ${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_DEBUG "-DBZ_DEBUG -g")

# -DNDEBUG?
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -Ofast -march=native")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -O3 -mllvm -vectorize -ffast-math -march=native") 
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "-Wfatal-errors -fast -axAVX")
endif()

# RELEASE mode by default
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

# from now on the 
find_package(OpenMP)

find_package(Boost COMPONENTS thread date_time system)
set(BOOST_THREAD_LIBS "${Boost_LIBRARIES}")

add_subdirectory(kahan_sum)
add_subdirectory(gnuplot-iostream_1d)
add_subdirectory(gnuplot-iostream_2d)
add_subdirectory(var_sign_2d)
add_subdirectory(harmosc)
add_subdirectory(todo_3d)
add_subdirectory(bombel)
add_subdirectory(concurrent_1d)
add_subdirectory(shallow_water)
add_subdirectory(rotating_cone)
add_subdirectory(particles)
add_subdirectory(lcph_kessler)
add_subdirectory(bcbf)
