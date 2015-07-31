# this is a compile-only test, nothing us run during "make test"
add_executable(mltpl_transl_unit test_mltpl_transl_unit.cpp u1.cpp u2.cpp)
target_link_libraries(mltpl_transl_unit ${libmpdataxx_LIBRARIES})                       
target_include_directories(mltpl_transl_unit PUBLIC ${libmpdataxx_INCLUDE_DIRS})
