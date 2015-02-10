#libmpdataxx_add_test(convergence_1d)
#add_test(
#  NAME convergence_1d_diff 
#  COMMAND bash -c "
#    for i in ${CMAKE_CURRENT_SOURCE_DIR}/*.txt.gz; do
#      echo AQQ;
#      ls -la \$i;
#      ls -la ${CMAKE_CURRENT_BINARY_DIR}/`basename \${i/.gz/}`;
#      echo BQQ;
#      zcmp \$i ${CMAKE_CURRENT_BINARY_DIR}/`basename \${i/.gz/}` || exit 1; 
#    done;
#  "
#)
add_test(
  NAME convergence_1d_plot
  COMMAND bash -c "
    for i in ${CMAKE_CURRENT_SOURCE_DIR}/*.txt.gz; do
      python ../../../../tests/tutorial/2_convergence_1d/plot.py ${CMAKE_CURRENT_BINARY_DIR}/`basename \${i/.gz/}`;
    done;
  "
)
