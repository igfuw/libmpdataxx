set(tmp test_concurrent_1d)
add_executable(${tmp} ${tmp}.cpp)

set_target_properties(${tmp} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  get_target_property(flags ${tmp} LINK_FLAGS)
  if (${flags} STREQUAL "flags-NOTFOUND") 
    set(flags "")
  endif()
  set_target_properties(${tmp} PROPERTIES LINK_FLAGS "${flags} ${OpenMP_CXX_FLAGS}")

  get_target_property(flags ${tmp} COMPILE_FLAGS)
  if (${flags} STREQUAL "flags-NOTFOUND")
    set(flags "")
  endif()
  set_target_properties(${tmp} PROPERTIES COMPILE_FLAGS "${flags} ${OpenMP_CXX_FLAGS}")
endif()

# Boost.Thread
find_package(Boost COMPONENTS thread date_time system REQUIRED)
target_link_libraries(${tmp} ${Boost_LIBRARIES})
find_package(Threads REQUIRED)
target_link_libraries(${tmp} ${CMAKE_THREAD_LIBS_INIT})

add_test(${tmp} ${tmp})
