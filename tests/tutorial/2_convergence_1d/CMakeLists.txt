#libmpdataxx_add_test(convergence_1d)
add_test(
  NAME convergence_1d_diff 
  COMMAND bash -c "
    for i in ${CMAKE_CURRENT_SOURCE_DIR}/*.txt.gz; 
      do zcmp \$i ${CMAKE_CURRENT_BINARY_DIR}/`basename \${i/.gz/}` || exit 1; 
    done;
  "
)
#add_test(NAME convergence_1d_plot) #TODO
