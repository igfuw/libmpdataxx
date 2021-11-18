libmpdataxx_add_test_gi(boussinesq_2d false) # don't compare results point-by-point with reference, because compiler optimization (-Ofast) affects them a little
if(NOT USE_MPI)
  add_test(
    NAME boussinesq_stats_diff 
    COMMAND zdiff ${CMAKE_CURRENT_SOURCE_DIR}/refdata/stats.txt.gz stats.txt
  )
endif()

set_property(TEST init_boussinesq_2d calc_boussinesq_2d clean_boussinesq_2d PROPERTY LABELS SlowWithMpi) 
