set(tmp test_todo_3d)
add_executable(${tmp} ${tmp}.cpp)
set_target_properties(${tmp} PROPERTIES LINK_FLAGS_DEBUG "-lblitz")
add_test(${tmp} ${tmp})
