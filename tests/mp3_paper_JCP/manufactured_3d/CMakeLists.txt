libmpdataxx_add_test(manufactured_3d)

add_test(manufactured_3d_stats_and_plots bash -c "
  python  ${CMAKE_CURRENT_SOURCE_DIR}/stats_and_plots.py out*;
")

add_test(manufactured_3d_conv_diff bash -c "
  diff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/manufactured_3d_conv.txt ${CMAKE_CURRENT_BINARY_DIR}/manufactured_3d_conv.txt;
")
