libmpdataxx_add_test(var_dt_test)
set_property(TEST var_dt_test PROPERTY LABELS NotWithMpi) # the test is super slow with MPI, make it easy to disable it
