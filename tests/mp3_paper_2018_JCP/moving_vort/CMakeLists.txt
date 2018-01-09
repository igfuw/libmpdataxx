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

add_test(moving_vort_diff bash -c "
  for dir in out_Mp3_192; do 
    echo   'comparing const.h5'                                                                                           &&
    h5diff --relative=2e-7 -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5               &&
    echo   'comparing temp.xmf'                                                                                           &&
    diff    $dir/temp.xmf                 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                               &&
    echo   'comparing timestep0000000000.h5'                                                                              &&
    h5diff --relative=2e-7 -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  &&
    echo   'comparing timestep0000000000.xmf'                                                                             &&
    diff    $dir/timestep0000000000.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                 &&
    echo   'comparing timestep0000047421.h5'                                                                              &&
    h5diff --relative=1e-5 -v $dir/timestep0000047421.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000047421.h5  &&
    echo   'comparing timestep0000047421.xmf'                                                                             &&
    diff    $dir/timestep0000047421.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000047421.xmf || exit 1;
  done 
")
