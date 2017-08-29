libmpdataxx_add_test(manufactured_3d)
add_test(manufactured__3d_stats_and_plots bash -c "
  python  ${CMAKE_CURRENT_SOURCE_DIR}/stats_and_plots.py mp3*;
") 
