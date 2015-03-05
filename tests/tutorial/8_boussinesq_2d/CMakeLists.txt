libmpdataxx_add_test_gi(boussinesq_2d)
add_test(
  NAME boussinesq_stats_diff 
  COMMAND zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats.txt.gz stats.txt
)
