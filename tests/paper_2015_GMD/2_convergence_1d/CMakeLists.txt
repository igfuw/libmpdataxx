libmpdataxx_add_test(convergence_1d)
if(NOT USE_MPI)
  add_test(
    NAME convergence_1d_diff 
    COMMAND bash -c "
      for i in ${CMAKE_CURRENT_SOURCE_DIR}/refdata/*.txt.gz; do
        echo \$i ;
        echo  compare to ; 
        echo ${CMAKE_CURRENT_BINARY_DIR}/`basename \${i/.gz/}`;
        zdiff \$i ${CMAKE_CURRENT_BINARY_DIR}/`basename \${i/.gz/}` || exit 1; 
      done;
    "
  )
  add_test(
    NAME convergence_1d_plot
    COMMAND bash -c "
      for i in ${CMAKE_CURRENT_SOURCE_DIR}/refdata/*.txt.gz; do
        python3 ${CMAKE_CURRENT_SOURCE_DIR}/plot.py ${CMAKE_CURRENT_BINARY_DIR}/`basename \${i/.gz/}`;
      done;
    "
  )
endif()
set_property(TEST convergence_1d PROPERTY LABELS SlowWithMpi) # the test is super slow with MPI, make it easy to disable it
