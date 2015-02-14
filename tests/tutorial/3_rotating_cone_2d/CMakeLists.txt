libmpdataxx_add_test_gi(rotating_cone_2d)
add_test(
  NAME rotating_cone_stats_diff 
  COMMAND bash -c "for i in stats*.txt; do zcmp ${CMAKE_CURRENT_SOURCE_DIR}/$i.gz $i || exit 1; done"
)
