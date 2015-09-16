libmpdataxx_add_test(revolving_sphere_3d)

add_test(revolving_sphere_3d_diff bash -c "
  for dir in basic fct iga iga_fct upwind; do 
    echo   'comparing const.h5'                                                                                        &&
    h5diff -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/const.h5                                 &&
    echo   'comparing temp.xmf'                                                                                        &&
    diff    $dir/temp.xmf               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                              &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5               &&
    echo   'comparing timestep0000000000.xmf'                                                                          &&
    diff    $dir/timestep0000000000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                &&
    echo   'comparing timestep0000000556.h5'                                                                           &&
    h5diff --delta=1e-16 -v $dir/timestep0000000556.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000556.h5 &&
    echo   'comparing timestep0000000556.xmf'                                                                          &&
    diff    $dir/timestep0000000556.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000556.xmf                || exit 1;
  done 
")

add_test(revolving_sphere_3d_stats_diff bash -c "
  for i in basic fct iga iga_fct upwind; do 
    echo ${CMAKE_CURRENT_SOURCE_DIR}/stats_$i.txt.gz;
    echo ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt;
    zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats_$i.txt.gz ${CMAKE_CURRENT_BINARY_DIR}/stats_$i.txt || exit 1; 
  done
") 

if(NOT $ENV{DISPLAY} STREQUAL "")
  add_test(revolving_sphere_3d_plot bash -c "
    for i in basic fct iga iga_fct upwind; do 
      pvpython ${CMAKE_CURRENT_SOURCE_DIR}/plot.py ${CMAKE_CURRENT_BINARY_DIR}/$i/timestep0000000556.xmf;
    done
  ")
else()
  message(STATUS "Skipping Paraview-based test (DISPLAY not set)")
endif()
#TODO: see https://public.kitware.com/Bug/view.php?id=13825
#set_tests_properties(revolving_sphere_3d_plot PROPERTIES SKIP_RETURN_CODE 44) 
