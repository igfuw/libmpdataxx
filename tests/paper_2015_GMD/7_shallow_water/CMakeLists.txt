libmpdataxx_add_test_hdf(shallow_water_1d)
add_test(shallow_water_1d_diff bash -c "
  for dir in 1d_fct_abs 1d_fct_iga; do 
    echo   'comparing const.h5'                                                                                        &&
    h5diff -v $dir/const.h5                            ${CMAKE_CURRENT_SOURCE_DIR}/refdata/const_1d.h5                 &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff -v $dir/timestep0000000000.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  &&
    echo   'comparing timestep0000000556.h5'                                                                           &&
    h5diff --delta=1e-18 -v $dir/timestep0000000300.h5 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000300.h5  || exit 1;
  done 
")
add_test(shallow_water_1d_plot bash -c "
  for dir in 1d_fct_abs 1d_fct_iga; do
    python  ${CMAKE_CURRENT_SOURCE_DIR}/plot_1d.py $dir;
  done
") 
add_test(shallow_water_1d_stats bash -c "
for dir in 1d_fct_abs 1d_fct_iga; do
    python  ${CMAKE_CURRENT_SOURCE_DIR}/errors_1d.py $dir                           &&
    zdiff $dir'_stats.txt' ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir'_stats.txt.gz'  || exit 1;
  done
") 
libmpdataxx_add_test_hdf(shallow_water_2d)
add_test(shallow_water_2d_diff bash -c "
  for dir in 2d_fct_abs 2d_fct_iga; do 
    echo   'comparing const.h5'                                                                                        &&
    h5diff -v $dir/const.h5                            ${CMAKE_CURRENT_SOURCE_DIR}/refdata/const_2d.h5                 &&
    echo   'comparing temp.xmf'                                                                                        &&
    diff    $dir/temp.xmf                              ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf               &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff -v $dir/timestep0000000000.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5  &&
    echo   'comparing timestep0000000000.xmf'                                                                          &&
    diff    $dir/timestep0000000000.xmf                ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf &&
    echo   'comparing timestep0000000556.h5'                                                                           &&
    h5diff --delta=2e-5 -v $dir/timestep0000000300.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000300.h5  &&
    echo   'comparing timestep0000000300.xmf'                                                                          &&
    diff    $dir/timestep0000000300.xmf                ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000300.xmf || exit 1;
  done 
")
 add_test(shallow_water_2d_plot bash -c "
  for dir in 2d_fct_abs 2d_fct_iga; do
    python  ${CMAKE_CURRENT_SOURCE_DIR}/plot_2d.py $dir;
  done
") 
add_test(shallow_water_3d_plot bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/plot_2d_3d.py;
") 
add_test(shallow_water_2d_stats bash -c "
  for dir in 2d_fct_abs 2d_fct_iga; do
    python  ${CMAKE_CURRENT_SOURCE_DIR}/errors_2d.py $dir                           && 
    zdiff $dir'_stats.txt' ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir'_stats.txt.gz'  || exit 1;
  done
") 
