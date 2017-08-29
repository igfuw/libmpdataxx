libmpdataxx_add_test(moving_vort)
add_test(moving_vort_stats_and_plots bash -c "
  python  ${CMAKE_CURRENT_SOURCE_DIR}/stats_and_plots.py out*;
")
if(FULL_SIM STREQUAL "true")
  add_test(moving_vort_conv_diff bash -c "
    diff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/moving_vort_conv.txt ${CMAKE_CURRENT_BINARY_DIR}/moving_vort_conv.txt;
  ")
else()
  add_test(moving_vort_conv_diff bash -c "
    diff <(head -n 7 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/moving_vort_conv.txt) ${CMAKE_CURRENT_BINARY_DIR}/moving_vort_conv.txt;
  ")
endif()
