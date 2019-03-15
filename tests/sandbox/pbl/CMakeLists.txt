libmpdataxx_add_test(pbl_iles)
libmpdataxx_add_test(pbl_smg)

add_test(pbl_iles_profiles bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_iles
")
add_test(pbl_smg_profiles bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_smg
")

add_test(pbl_iles_budget bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_iles
")
add_test(pbl_smg_budget bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_smg
")

libmpdataxx_add_test_travis(pbl_smg_travis)
libmpdataxx_add_test_travis(pbl_iles_travis)

add_test(pbl_iles_travis_profiles bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_iles_travis
")
add_test(pbl_smg_travis_profiles bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_smg_travis
")

add_test(pbl_iles_travis_prof_diff
        zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/profiles_pbl_iles_travis.txt.gz profiles_pbl_iles_travis.txt)
add_test(pbl_smg_travis_prof_diff
        zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/profiles_pbl_smg_travis.txt.gz profiles_pbl_smg_travis.txt)

add_test(pbl_iles_travis_budget bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_iles_travis
")
add_test(pbl_smg_travis_budget bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_smg_travis
")

add_test(pbl_iles_travis_budget_diff
        zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/budget_pbl_iles_travis.txt.gz budget_pbl_iles_travis.txt)
add_test(pbl_smg_travis_budget_diff
        zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/budget_pbl_smg_travis.txt.gz budget_pbl_smg_travis.txt)

add_test(pbl_iles_travis_diff bash -c "
    dir=out_pbl_iles_travis
    echo   'comparing const.h5'                                                                                        &&
    h5diff -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5                            &&
    echo   'comparing temp.xmf'                                                                                        &&
    diff    $dir/temp.xmf               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                              &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5               &&
    echo   'comparing timestep0000000000.xmf'                                                                          &&
    diff    $dir/timestep0000000000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                &&
    echo   'comparing timestep0000000600.h5'                                                                           &&
    h5diff --delta=1e-5 -v $dir/timestep0000000600.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.h5  &&
    echo   'comparing timestep0000000600.xmf'                                                                          &&
    diff    $dir/timestep0000000600.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.xmf                || exit 1;
")

add_test(pbl_smg_travis_diff bash -c "
    dir=out_pbl_smg_travis
    echo   'comparing const.h5'                                                                                        &&
    h5diff -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5                            &&
    echo   'comparing temp.xmf'                                                                                        &&
    diff    $dir/temp.xmf               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                              &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5               &&
    echo   'comparing timestep0000000000.xmf'                                                                          &&
    diff    $dir/timestep0000000000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                &&
    echo   'comparing timestep0000000600.h5'                                                                           &&
    h5diff --delta=1e-5 -v $dir/timestep0000000600.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.h5  &&
    echo   'comparing timestep0000000600.xmf'                                                                          &&
    diff    $dir/timestep0000000600.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.xmf                || exit 1;
")
