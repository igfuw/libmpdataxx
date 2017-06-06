libmpdataxx_add_test(moving)
set(comp_body "
    echo   'comparing const.h5'                                                                                         &&
    h5diff -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/const.h5                                  &&
    echo   'comparing temp.xmf'                                                                                         &&
    diff    $dir/temp.xmf                 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                             &&
    echo   'comparing timestep0000000000.h5'                                                                            &&
    h5diff -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5                &&
    echo   'comparing timestep0000000000.xmf'                                                                           &&
    diff    $dir/timestep0000000000.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf               &&
    echo   'comparing timestep0000006004.h5'                                                                            &&
    h5diff --delta=1e-18 -v $dir/timestep0000006004.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000006004.h5  &&
    echo   'comparing timestep0000006004.xmf'                                                                           &&
    diff    $dir/timestep0000006004.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000006004.xmf               &&
    echo   'comparing timestep0000012008.h5'                                                                            &&
    h5diff --delta=1e-18 -v $dir/timestep0000012008.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000012008.h5  &&
    echo   'comparing timestep0000012008.xmf'                                                                           &&
    diff    $dir/timestep0000012008.xmf   ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000012008.xmf || exit 1;
")

add_test(moving_diff bash -c "
  for dir in nug_i2_96 nug_iga_fct_i2_96; do 
  ${comp_body}
  done 
")

if(NOT USE_MPI)
  add_test(moving_stats_diff bash -c "
    for i in nug_i2 nug_iga_fct_i2; do 
      echo  ${CMAKE_CURRENT_SOURCE_DIR}/stats_$i.txt.gz;
      echo  ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt;
      zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats_$i.txt.gz ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt || exit 1; 
    done
  ") 
endif()

libmpdataxx_add_test(moving_tot)

add_test(moving_tot_diff bash -c "
  for dir in nug_tot_i3_96 nug_iga_tot_fct_i2_96; do 
  ${comp_body}
  done 
")

if(NOT USE_MPI)
  add_test(moving_tot_stats_diff bash -c "
    for i in nug_tot_i3 nug_iga_tot_fct_i2; do 
      echo  ${CMAKE_CURRENT_SOURCE_DIR}/stats_$i.txt.gz;
      echo  ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt;
      zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats_$i.txt.gz ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt || exit 1; 
    done
  ") 
endif()
