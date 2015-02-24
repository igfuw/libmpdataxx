libmpdataxx_add_test_gi(coupled_harmosc)
add_test(
  NAME harmosc_stats_diff 
  COMMAND zdiff ${CMAKE_CURRENT_SOURCE_DIR}/stats.txt.gz stats.txt
)
