libmpdataxx_add_test(stationary)

add_test(stationary_diff bash -c "
  for dir in nug_i2_96 nug_iga_fct_i2_96 nug_tot_i3_96 nug_iga_tot_fct_i2_96; do 
    echo   'comparing const.h5'                                                                                         &&
    h5diff -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/const.h5                                  &&
    echo   'comparing temp.xmf'                                                                                         &&
    diff    $dir/temp.xmf                 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                             &&
    echo   'comparing timestep0000000000.h5'                                                                            &&
    h5diff -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5                &&
    echo   'comparing timestep0000000000.xmf'                                                                           &&
    diff    $dir/timestep0000000000.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf               &&
    echo   'comparing timestep0000001492.h5'                                                                            &&
    h5diff --delta=1e-18 -v $dir/timestep0000001492.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000001492.h5  &&
    echo   'comparing timestep0000001492.xmf'                                                                           &&
    diff    $dir/timestep0000001492.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000001492.xmf               &&
    echo   'comparing timestep0000002983.h5'                                                                            &&
    h5diff --delta=1e-18 -v $dir/timestep0000002983.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000002983.h5  &&
    echo   'comparing timestep0000002983.xmf'                                                                           &&
    diff    $dir/timestep0000002983.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000002983.xmf || exit 1;
  done 
")

add_test(stationary_stats_diff bash -c "
  for i in nug_i2 nug_iga_fct_i2 nug_tot_i3 nug_iga_tot_fct_i2; do 
    echo  ${CMAKE_CURRENT_SOURCE_DIR}/stats_$i.txt.gz;
    echo  ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt;
    zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats_$i.txt.gz ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt || exit 1; 
  done
") 
