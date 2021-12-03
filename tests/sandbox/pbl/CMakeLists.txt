libmpdataxx_add_test(pbl_iles)
libmpdataxx_add_test(pbl_smg)

add_test(pbl_iles_profiles bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_iles
")
add_test(pbl_smg_profiles bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_smg
")

add_test(pbl_iles_budget bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_iles
")
add_test(pbl_smg_budget bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_smg
")

libmpdataxx_add_test(pbl_smg_short)
libmpdataxx_add_test(pbl_iles_short)

add_test(pbl_iles_short_profiles bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_iles_short
")
add_test(pbl_smg_short_profiles bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py out_pbl_smg_short
")

add_test(pbl_iles_short_prof_diff
        zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/profiles_pbl_iles_short.txt.gz profiles_pbl_iles_short.txt)
add_test(pbl_smg_short_prof_diff
        zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/profiles_pbl_smg_short.txt.gz profiles_pbl_smg_short.txt)

add_test(pbl_iles_short_budget bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_iles_short
")
add_test(pbl_smg_short_budget bash -c "
    python3  ${CMAKE_CURRENT_SOURCE_DIR}/budget.py out_pbl_smg_short
")

add_test(pbl_iles_short_budget_diff
        h5diff -v2 --relative=1e-2 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/budget_pbl_iles_short.h5 budget_pbl_iles_short.h5)
add_test(pbl_smg_short_budget_diff
        h5diff -v2 --relative=1e-2 ${CMAKE_CURRENT_SOURCE_DIR}/refdata/budget_pbl_smg_short.h5 budget_pbl_smg_short.h5)

add_test(pbl_iles_short_diff bash -c "
    dir=out_pbl_iles_short
    echo   'comparing const.h5'                                                                                        &&
    h5diff --use-system-epsilon -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5                            &&
    echo   'comparing temp.xmf'                                                                                        &&
    diff    $dir/temp.xmf               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                              &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff --use-system-epsilon -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5               &&
    echo   'comparing timestep0000000000.xmf'                                                                          &&
    diff    $dir/timestep0000000000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                &&
    echo   'comparing timestep0000000600.h5'                                                                           &&
    h5diff --delta=1e-5 -v $dir/timestep0000000600.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.h5  &&
    echo   'comparing timestep0000000600.xmf'                                                                          &&
    diff    $dir/timestep0000000600.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.xmf                || exit 1;
")

add_test(pbl_smg_short_diff bash -c "
    dir=out_pbl_smg_short
    echo   'comparing const.h5'                                                                                        &&
    h5diff --use-system-epsilon -v $dir/const.h5               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/const.h5                            &&
    echo   'comparing temp.xmf'                                                                                        &&
    diff    $dir/temp.xmf               ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/temp.xmf                              &&
    echo   'comparing timestep0000000000.h5'                                                                           &&
    h5diff --use-system-epsilon -v $dir/timestep0000000000.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.h5               &&
    echo   'comparing timestep0000000000.xmf'                                                                          &&
    diff    $dir/timestep0000000000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000000.xmf                &&
    echo   'comparing timestep0000000600.h5'                                                                           &&
    h5diff --delta=1e-5 -v $dir/timestep0000000600.h5  ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.h5  &&
    echo   'comparing timestep0000000600.xmf'                                                                          &&
    diff    $dir/timestep0000000600.xmf ${CMAKE_CURRENT_SOURCE_DIR}/refdata/$dir/timestep0000000600.xmf                || exit 1;
")

# test labels
#set_property(TEST pbl_iles       pbl_smg       pbl_iles_budget       pbl_smg_budget       pbl_iles_profiles       pbl_smg_profiles       PROPERTY LABELS NotInCI)
set_property(TEST pbl_iles       pbl_smg       pbl_iles_budget       pbl_smg_budget       pbl_iles_profiles       pbl_smg_profiles       PROPERTY LABELS SlowWithMpi)
set_property(TEST pbl_iles_short pbl_smg_short pbl_iles_short_budget pbl_smg_short_budget pbl_iles_short_profiles pbl_smg_short_profiles PROPERTY LABELS SlowWithMpi ShorterVersion)
set_property(TEST pbl_iles_short_diff pbl_smg_short_diff pbl_iles_short_budget_diff pbl_smg_short_budget_diff pbl_iles_short_prof_diff pbl_smg_short_prof_diff PROPERTY LABELS SlowWithMpi ShorterVersion)
