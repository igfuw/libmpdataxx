libmpdataxx_add_test_gi(coupled_harmosc)
if(NOT USE_MPI)
  add_test(
    NAME harmosc_stats_diff 
    COMMAND zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats.txt.gz stats.txt
  )
endif()
