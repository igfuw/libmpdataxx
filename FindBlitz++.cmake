find_path(BLITZ_INCLUDE_DIR blitz/array.h)
find_library(BLITZ_LIBRARY NAMES blitz)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Blitz DEFAULT_MSG BLITZ_LIBRARY BLITZ_INCLUDE_DIR)

include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_INCLUDES ${BLITZ_INCLUDE_DIR})
check_include_file_cxx("blitz/array.h" BLITZ_FOUND)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "Blitz++ includes not found - please install Blitz++ or point CMake to it.")
endif()
check_cxx_source_compiles("#include <blitz/array.h>\nint main(){blitz::Array<float,1> a(2); a = blitz::safeToReturn(a+a);}" BLITZ_VER_AT_LEAST_0_10)
if (NOT BLITZ_VER_AT_LEAST_0_10)
  message(FATAL_ERROR "Blitz++ ver >= 0.10 requirement not met - please update Blitz++.")
endif()

