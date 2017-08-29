libmpdataxx_add_test(reversing_deform)
add_test(reversing_deform_stats_and_plots bash -c "
  python  ${CMAKE_CURRENT_SOURCE_DIR}/stats_and_plots.py out*;
") 
if(FULL_SIM STREQUAL "true")
  add_test(reversing_deform_conv_diff bash -c "
    diff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/reversing_deform_conv.txt ${CMAKE_CURRENT_BINARY_DIR}/reversing_deform_conv.txt;
  ")
else()
  add_test(reversing_deform_conv_diff bash -c "
    diff <(head -n 6 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/reversing_deform_conv.txt) ${CMAKE_CURRENT_BINARY_DIR}/reversing_deform_conv.txt;
  ")
endif()
