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
    diff <(head -n 7 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/reversing_deform_conv.txt) ${CMAKE_CURRENT_BINARY_DIR}/reversing_deform_conv.txt;
  ")
endif()

add_test(reversing_deform_mixing_diff bash -c "
  diff <(head -n 7 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/reversing_deform_mixing.txt) ${CMAKE_CURRENT_BINARY_DIR}/reversing_deform_mixing.txt;
")

add_test(reversing_deform_diff_120 bash -c "
  for dir in out_Mg3No_120; do 
    echo   'comparing const.h5'                                                                                           &&
    h5diff --relative=2e-7 -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5               &&
    echo   'comparing temp.xmf'                                                                                           &&
    diff    $dir/temp.xmf                 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                               &&
    echo   'comparing timestep0000000000.h5'                                                                              &&
    h5diff --relative=2e-7 -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  &&
    echo   'comparing timestep0000000000.xmf'                                                                             &&
    diff    $dir/timestep0000000000.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                 &&
    echo   'comparing timestep0000000454.h5'                                                                              &&
    h5diff --relative=2e-3 -v $dir/timestep0000000454.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000454.h5  &&
    echo   'comparing timestep0000000454.xmf'                                                                             &&
    diff    $dir/timestep0000000454.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000454.xmf                 &&
    echo   'comparing timestep0000000908.h5'                                                                              &&
    h5diff --relative=2e-3 -v $dir/timestep0000000908.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000908.h5  &&
    echo   'comparing timestep0000000908.xmf'                                                                             &&
    diff    $dir/timestep0000000908.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000908.xmf || exit 1;
  done 
")

add_test(reversing_deform_diff_240 bash -c "
  for dir in out_Mp3_240 out_Mp3cc_240 out_Mg2No_240 out_Mg3No_240; do 
    echo   'comparing const.h5'                                                                                           &&
    h5diff --relative=2e-7 -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5               &&
    echo   'comparing temp.xmf'                                                                                           &&
    diff    $dir/temp.xmf                 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                               &&
    echo   'comparing timestep0000000000.h5'                                                                              &&
    h5diff --relative=2e-7 -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  &&
    echo   'comparing timestep0000000000.xmf'                                                                             &&
    diff    $dir/timestep0000000000.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                 &&
    echo   'comparing timestep0000000908.h5'                                                                              &&
    h5diff --relative=2e-3 -v $dir/timestep0000000908.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000908.h5  &&
    echo   'comparing timestep0000000908.xmf'                                                                             &&
    diff    $dir/timestep0000000908.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000908.xmf                 &&
    echo   'comparing timestep0000001816.h5'                                                                              &&
    h5diff --relative=2e-3 -v $dir/timestep0000001816.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000001816.h5  &&
    echo   'comparing timestep0000001816.xmf'                                                                             &&
    diff    $dir/timestep0000001816.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000001816.xmf || exit 1;
  done 
")
