set(tmp test_kahan_sum)
add_executable(${tmp} ${tmp}.cpp)
set_target_properties(${tmp} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")
add_test(${tmp} ${tmp})
