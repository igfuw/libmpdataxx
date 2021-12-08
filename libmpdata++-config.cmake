if(APPLE)
  # needed for the XCode clang to be identified as AppleClang and not Clang
  cmake_minimum_required(VERSION 3.0) 
else()
  # needed for the OpenMP test to work in C++-only project 
  # (see http://public.kitware.com/Bug/view.php?id=11910)
  cmake_minimum_required(VERSION 2.8.8) 
endif()

cmake_minimum_required(VERSION 3.14) # for setting output variable of try_compile 

# the policies we care about:
# - CMP0025 - make CMake distinguis between Apple and LLVM clang
# - CMP0042 - make CMake use RPATHs on OSX
if(CMAKE_VERSION VERSION_GREATER 2.9)
  cmake_policy(VERSION 3.0)
endif()

############################################################################################
# the following variables will be set:
set(libmpdataxx_FOUND False)
set(libmpdataxx_INCLUDE_DIRS "")
set(libmpdataxx_LIBRARIES "")
set(libmpdataxx_CXX_FLAGS_DEBUG "")
set(libmpdataxx_CXX_FLAGS_RELEASE "")
set(libmpdataxx_MPIRUN "")

############################################################################################
# libmpdata++ headers for non-default install location (i.e. for make DESTDIR=<dir> install)
set(libmpdataxx_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../include/")


############################################################################################
# debug mode compiler flags
set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -std=c++14 -DBZ_DEBUG -g -Wno-enum-compare") #TODO: -Og if compiler supports it?


############################################################################################
# release mode compiler flags
if(
  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR 
  CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
  CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" 
)
  set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -std=c++14 -DNDEBUG -Ofast -march=native -Wno-enum-compare")

  # preventing Kahan summation from being optimised out
  if (
    (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6) OR
    (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1) #TODO: never actually checked!
  )
    set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -fno-vectorize") 
  endif()
endif()


if(
  CMAKE_CXX_COMPILER_ID STREQUAL "Intel"
)
  # flags taken from -fast but without -static
  set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -std=gnu++14 -DNDEBUG -xHOST -O3 -ipo -no-prec-div -fp-model fast=2")
endif()


############################################################################################
# C++14
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-std=c++14")
check_cxx_source_compiles("
  #include <type_traits>
  auto f() { return 1;}
  template <bool a, class b> using ei=std::enable_if<a,b>; 
  struct a {a(int){}};struct b:a {using a::a;};  
  int main(){b i(1);}
" CXX14_SUPPORTED)
if (NOT CXX14_SUPPORTED)
  message(FATAL_ERROR "C++14 compatibility test failed - please update your compiler or point CMake to another one with -DCMAKE_CXX_COMPILER=...")
endif()
unset(CMAKE_REQUIRED_FLAGS)


############################################################################################
# Blitz++
find_package(PkgConfig)
pkg_check_modules(BLITZ QUIET blitz>=0.10)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "Blitz++ library not found or Blitz++ version requirement not met (>=0.10)

* To insall Blitz++, please try:
*   Debian/Ubuntu: sudo apt-get install libblitz0-dev
*   Fedora: sudo yum install blitz-devel
*   Homebrew: brew install blitz
  ")
else()
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${BLITZ_INCLUDE_DIRS}")
  set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${BLITZ_LIBRARIES}")
endif()

############################################################################################
# OpenMP
find_package(OpenMP QUIET)
if(OPENMP_FOUND)
  set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
  set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
else()
  message(STATUS "OpenMP not supported by the compiler.

* Programs using libmpdata++'s OpenMP concurrency will compile but run in serial mode
  ")
endif()


############################################################################################
# multi-threading
# find_package(ThreadsCXX) <- this requires C language to be enabled
# TODO: better solution! 
# TODO: not needed for serial-only programs!
# TODO: -fopenmp implies -pthread on gcc
set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -pthread")
set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -pthread")


############################################################################################
# MPI - detecting if the C++ compiler is actually an MPI wrapper
set(msg "Detecting if the compiler is an MPI wrapper...")
message(STATUS "${msg}")
execute_process(COMMAND ${CMAKE_CXX_COMPILER} "-show" RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_QUIET)
if (status EQUAL 0 AND output MATCHES "mpi") 
  set(USE_MPI TRUE)
  set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -DUSE_MPI")
  set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -DUSE_MPI")
  set(libmpdataxx_MPIRUN ${CMAKE_CXX_COMPILER})
  string(REPLACE "mpic++" "mpirun" libmpdataxx_MPIRUN ${libmpdataxx_MPIRUN})
  string(REPLACE "mpicxx" "mpirun" libmpdataxx_MPIRUN ${libmpdataxx_MPIRUN})
  string(REPLACE "mpiXX"  "mpirun" libmpdataxx_MPIRUN ${libmpdataxx_MPIRUN})
else()
  set(USE_MPI FALSE)
endif()
message(STATUS "${msg} - ${USE_MPI}")
unset(msg)
unset(status)
unset(output)

############################################################################################
# Boost libraries
set(Boost_DETAILED_FAILURE_MSG ON)
set(req_comp thread date_time system iostreams timer filesystem)
if(USE_MPI)
  list(APPEND req_comp mpi)
  list(APPEND req_comp serialization)
  #set(Boost_VERSION 1.59.0)
else()
  # Boost libraries v>=1.55.0, because boost/predef was added then
  #set(Boost_VERSION 1.55.0)
endif()
find_package(Boost COMPONENTS ${req_comp})
if(Boost_FOUND)
  set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${Boost_LIBRARIES}")
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${Boost_INCLUDE_DIRS}")
  if(
    (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") 
    AND (Boost_MINOR_VERSION EQUAL 55) 
  )
    # add a definition -DBOOST_HAS_INT128=1 to clang calls on linux to avoid errors with boost.atomic (https://svn.boost.org/trac/boost/ticket/9610)
    set(libmpdataxx_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -DBOOST_HAS_INT128=1")
    set(libmpdataxx_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE} -DBOOST_HAS_INT128=1")
  endif()
else()
  #TODO: check separately for optional and mandatory components
  message(FATAL_ERROR "Boost (or some of its components) not found.

* Programs based on libmpdata++ will not compile. 
* To insall Boost, please try:
*   Debian/Ubuntu: sudo apt-get install libboost-all-dev
*   Fedora: sudo yum install boost-devel
  ")
endif()


############################################################################################
# HDF5 libraries
find_package(HDF5 COMPONENTS CXX HL)
if(HDF5_FOUND)
  if(NOT HDF5_CXX_LIBRARIES)
    message(FATAL_ERROR "HDF5 installation lacks C++ support.")
  endif()

  if(USE_MPI AND NOT HDF5_IS_PARALLEL)
    message(STATUS "MPI was enabled for libmpdata++ but not in HDF5.

* Programs using libmpdata++'s HDF5 output will not compile.
* To install MPI-enabled HDF5, please try:
*   Debian/Ubuntu: sudo apt-get install libhdf5-openmpi-dev
*   Fedora: sudo yum install hdf5-openmpi-devel
*   Homebrew: brew install hdf5 --with-cxx --with-mpi
    ") 
  endif()

  if(NOT USE_MPI AND HDF5_IS_PARALLEL)
    message(STATUS "MPI was enabled in HDF5 but not in libmpdata++.

* Programs using libmpdata++'s HDF5 output will not compile.
* To install serial HDF5, please try:
*   Debian/Ubuntu: sudo apt-get install libhdf5-serial-dev hdf5-tools  (TODO)
*   Fedora: sudo yum install hdf5-devel                                (TODO)
*   Homebrew: brew install hdf5 --with-cxx                             (TODO)
*
* To enable MPI in libmpdata++ point it to a MPI compiler 
* with -DCMAKE_CXX_COMPILER=<mpi_compiler>
  ")

  endif()

  if(USE_MPI AND HDF5_IS_PARALLEL)
    # detecting if HDF5-MPI is usable from C++
    # see http://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2015-June/008600.html
    # https://github.com/live-clones/hdf5/commit/cec2478e71d2358a2df32b3dbfeed8b0b51980bb
    set(msg "Checking if MPI-HDF5 is usable from C++...")
    set(pfx "HDF5/MPI/C++ check")
    message(STATUS ${msg})
    execute_process(COMMAND "mktemp" "-d" RESULT_VARIABLE status OUTPUT_VARIABLE tmpdir OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT status EQUAL 0)                                                       
      message(FATAL_ERROR "${pfx}: mkdtemp failed")                               
    endif()                                                                       
    file(WRITE "${tmpdir}/test.cpp" "                                              
      #include <boost/mpi/environment.hpp>
      #include <H5Cpp.h>
      #if !defined(H5_HAVE_PARALLEL)
      #  error H5_HAVE_PARALLEL not defined!
      #endif
      int main() 
      { 
        boost::mpi::environment e; 
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);                                
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);                  
        H5::H5File(\"test.h5\", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);                
        H5Pclose(plist_id); 
      }
    ")

#    message("${CMAKE_CXX_COMPILER}" "test.cpp" "-I${Boost_INCLUDE_DIRS}" "-I${HDF5_INCLUDE_DIRS}" ${HDF5_LIBRARIES} ${Boost_LIBRARIES})

#    execute_process(
#      COMMAND "${CMAKE_CXX_COMPILER}" "test.cpp" "-I${Boost_INCLUDE_DIRS}" "-I${HDF5_INCLUDE_DIRS}" ${HDF5_LIBRARIES} ${Boost_LIBRARIES}# the order of HDF/Boost matters here!
#      # Boost_LIBRARY_DIRS has to be in LD_RUN_PATH, tried to specify it through rpath/rpath-link but failed; at runtime it correctly finds libboost-mpi (direct dependency), but fails on boost-serialization (which is a dependency of boost-mpi)
#      #COMMAND "${CMAKE_CXX_COMPILER}" "test.cpp" "-I${Boost_INCLUDE_DIRS}" "-Wl,-rpath,${Boost_LIBRARY_DIRS},-rpath-link,${Boost_LIBRARY_DIRS}" ${HDF5_LIBRARIES} ${Boost_LIBRARIES}# the order of HDF/Boost matters here!
#      WORKING_DIRECTORY ${tmpdir} 
#      RESULT_VARIABLE status 
#      ERROR_VARIABLE error
#    )

    try_compile(status ${CMAKE_BINARY_DIR} "${tmpdir}/test.cpp" 
                OUTPUT_VARIABLE error
                CMAKE_FLAGS INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS},{HDF5_INCLUDE_DIRS}
                LINK_LIBRARIES ${Boost_LIBRARIES}
                LINK_LIBRARIES ${HDF5_LIBRARIES}
                COPY_FILE ${tmpdir}/test_hdf5_mpi
                COPY_FILE_ERROR copy_error
               )

    if (status EQUAL 0)                                                       
      message(FATAL_ERROR "${pfx}: compilation failed\n status: ${status}\n copy file error: ${copy_error}\n output: ${error}")                               
    endif()                                                                       
    message(STATUS "${msg} - compilation OK")

    execute_process(
      COMMAND "mpiexec" "-np" "1" "./test_hdf5_mpi" 
      WORKING_DIRECTORY ${tmpdir} 
      RESULT_VARIABLE status
      ERROR_VARIABLE error
    )
    if (NOT status EQUAL 0)                                                       
      message(FATAL_ERROR "${pfx}: execution failed\n ${error}
        likely you have to upgrade HDF5, see:
        - http://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2015-June/008600.html
        - https://github.com/live-clones/hdf5/commit/cec2478e71d2358a2df32b3dbfeed8b0b51980bb
      ")
    endif()                                                                       
    message(STATUS "${msg} - non-mpirun execution OK")

    # detecting if it runs under mpirun (missing libhwloc-plugins issue:
    # https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=790540
    # )
    execute_process(COMMAND ${libmpdataxx_MPIRUN} "-np" "2" "./test_hdf5_mpi" 
      WORKING_DIRECTORY ${tmpdir} 
      RESULT_VARIABLE status
      ERROR_VARIABLE error
    )
    if (NOT status EQUAL 0)                                                       
      message(FATAL_ERROR "TODO: ${status}\n ${error}")
    endif()
    message(STATUS "${msg} - mpirun execution OK")

    unset(status)
    unset(pfx)
    unset(msg)
  endif() 


  #
  set(libmpdataxx_LIBRARIES "${libmpdataxx_LIBRARIES};${HDF5_LIBRARIES}")
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${HDF5_INCLUDE_DIRS}")
else()
  message(STATUS "HDF5 not found. 

* Programs using libmpdata++'s HDF5 output will not compile.
* To install HDF5, please try:
*   Debian/Ubuntu: sudo apt-get install libhdf5-serial-dev hdf5-tools
*   Fedora: sudo yum install hdf5-devel
*   Homebrew: brew install hdf5 --with-cxx
  ")
endif()


############################################################################################
# gnuplot-iostream
find_path(GNUPLOT-IOSTREAM_INCLUDE_DIR PATH_SUFFIXES gnuplot-iostream/ NAMES gnuplot-iostream.h)
if(GNUPLOT-IOSTREAM_INCLUDE_DIR)
  set(libmpdataxx_INCLUDE_DIRS "${libmpdataxx_INCLUDE_DIRS};${GNUPLOT-IOSTREAM_INCLUDE_DIR}")
else()
  message(STATUS "gnuplot-iostream not found.

* Programs using libmpdata++'s gnuplot-iostream output will not compile.
* To install gnuplot-iostream, please try:
*   Debian/Ubuntu: sudo apt-get install libgnuplot-iostream-dev
*   manual: wget -O /usr/local/include/gnuplot-iostream.h https://raw.githubusercontent.com/dstahlke/gnuplot-iostream/master/gnuplot-iostream.h
  ")
endif()

find_program(GNUPLOT_FOUND NAMES gnuplot)
if(GNUPLOT_FOUND)
else()
  message(STATUS "gnuplot not found.

* Programs using libmpdata++'s gnuplot-iostream output will not run.
* To install gnuplot, please try:
*   Debian/Ubuntu: sudo apt-get install gnuplot
*   Fedora: sudo yum install gnuplot
*   Homebrew: brew install gnuplot
  ")
endif()


############################################################################################
list(REMOVE_DUPLICATES libmpdataxx_INCLUDE_DIRS)
list(REMOVE_ITEM libmpdataxx_INCLUDE_DIRS "")
set(libmpdataxx_FOUND TRUE)
