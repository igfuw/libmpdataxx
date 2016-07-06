libmpdataxx_add_test(pbl)
add_test(pbl_profiles bash -c "
    python  ${CMAKE_CURRENT_SOURCE_DIR}/profiles.py ${CMAKE_CURRENT_BINARY_DIR}/out_pbl/timestep0000001500.h5
")
