libmpdataxx_add_test(over_the_pole_2d)
add_test(over_the_pole_2d_diff bash -c "
  for dir in best default; do 
    echo   'comparing coord.h5'                                                                        &&
    h5diff $dir/coord.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/coord.h5               &&
    echo   'comparing temp.xmf'                                                                        &&
    cmp    $dir/temp.xmf               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf               &&
    echo   'comparing timestep0000000000.h5'                                                           &&
    h5diff $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  &&
    echo   'comparing timestep0000000000.xmf'                                                          &&
    cmp    $dir/timestep0000000000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf &&
    echo   'comparing timestep0000005120.h5'                                                           &&
    h5diff $dir/timestep0000005120.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000005120.h5  &&
    echo   'comparing timestep0000005120.xmf'                                                          &&
    cmp    $dir/timestep0000005120.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000005120.xmf;
  done 
")
