# copy of function defined in tests/paper_2015_GMD/CMakeLists.txt ...
function(libmpdataxx_add_test_gi test diff)
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${libmpdataxx_LIBRARIES})
  target_include_directories(${test} PUBLIC ${libmpdataxx_INCLUDE_DIRS})

  # this will not work with multiple threads (make -j) as it depends on the order of execution + $PPID might be different?
  add_test(NAME "init_${test}" COMMAND "bash" "-c" "> log-$PPID")
  if(USE_MPI)
    add_test(NAME "calc_${test}" COMMAND "bash" "-c" "GNUPLOT_IOSTREAM_CMD=\"cat >> log-`echo $PPID`\" ${libmpdataxx_MPIRUN} -np 3 ${CMAKE_CURRENT_BINARY_DIR}/${test}")
  else()
    add_test(NAME "calc_${test}" COMMAND "bash" "-c" "GNUPLOT_IOSTREAM_CMD=\"cat >> log-`echo $PPID`\" ${CMAKE_CURRENT_BINARY_DIR}/${test}")
  endif()
  if(NOT USE_MPI)
    add_test(NAME "plot_${test}" COMMAND "bash" "-c" "gnuplot log-$PPID")
    if(${diff})
      add_test(NAME "diff_${test}" COMMAND "bash" "-c" "zdiff log-$PPID ${CMAKE_CURRENT_SOURCE_DIR}/refdata/log-orig.gz && rm log-$PPID")
    else()
     add_test(NAME "clean_${test}" COMMAND "bash" "-c" "rm log-$PPID")
    endif()
  else()
    add_test(NAME "clean_${test}" COMMAND "bash" "-c" "rm log-$PPID")
  endif()
endfunction()

libmpdataxx_add_test_gi(inflow_2d false) # don't compare results point-by-point with reference, because compiler optimization (-Ofast) affects them a little

set_property(TEST init_inflow_2d calc_inflow_2d PROPERTY LABELS SlowWithMpi) 
