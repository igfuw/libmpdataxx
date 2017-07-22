libmpdataxx_add_test(pbl_iles)
libmpdataxx_add_test(pbl_smg)

libmpdataxx_add_test_travis(pbl_smg_travis)
libmpdataxx_add_test_travis(pbl_iles_travis)

#add_test(pbl_iles_profiles bash -c "
#    python  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py ${CMAKE_CURRENT_BINARY_DIR}/out_pbl_iles/timestep0000001500.h5
#")
#add_test(pbl_smg_profiles bash -c "
#    python  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py ${CMAKE_CURRENT_BINARY_DIR}/out_pbl_smg/timestep0000001500.h5
#")
