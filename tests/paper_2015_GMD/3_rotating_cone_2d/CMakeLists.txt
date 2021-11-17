libmpdataxx_add_test_gi(rotating_cone_2d false) # disable comparison of data sent to Gnuplot, because the data depends on compiler's optimization (result of blitz::rint in particular?)
if(NOT USE_MPI)
  add_test(
    NAME rotating_cone_stats_diff 
    COMMAND bash -c "for i in stats*.txt; do 
      zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$i.gz $i || exit 1; 
    done"
  )
endif()
