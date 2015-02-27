libmpdataxx_add_test(revolving_sphere_3d)
add_test(revolving_sphere_3d_diff bash -c "
  for dir in basic fct iga iga_fct upwind; do 
    echo   'comparing coord.h5'                                                                                        &&
    h5diff -v $dir/coord.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/coord.h5                            &&
    echo   'comparing temp.xmf'                                                                                        &&
    diff    $dir/temp.xmf               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                              &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5               &&
    echo   'comparing timestep0000000000.xmf'                                                                          &&
    diff    $dir/timestep0000000000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                &&
    echo   'comparing timestep0000000556.h5'                                                                           &&
    h5diff --delta=1e-17 -v $dir/timestep0000000556.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000556.h5 &&
    echo   'comparing timestep0000000556.xmf'                                                                          &&
    diff    $dir/timestep0000000556.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000556.xmf || exit 1;
  done 
")
add_test(revolving_sphere_3d_stats_diff bash -c "
  for i in basic fct iga iga_fct upwind; do 
    echo ${CMAKE_CURRENT_SOURCE_DIR}/stats_$i.txt.gz;
    echo ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt;
    zdiff ${CMAKE_CURRENT_SOURCE_DIR}/stats_$i.txt.gz ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt || exit 1; 
  done
") 

