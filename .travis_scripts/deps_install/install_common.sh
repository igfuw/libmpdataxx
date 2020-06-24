#!/usr/bin/env sh
set -e
############################################################################
## All the cached dependencies are installed in ${TRAVIS_BUILD_DIR}/deps/
#############################################################################
DEPS_DIR="${TRAVIS_BUILD_DIR}/deps"

# silence the gazillion warnings coming from blitz headers when using the osx clang
if [[ $TRAVIS_OS_NAME == 'osx' && $COMPILER == 'clang++' ]]; then export CXXFLAGS="-Wno-parentheses ${CXXFLAGS}"; fi

# get libclang-dev for headers
if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'clang++' && $MPI != 'none' ]]; then sudo $apt_get_install libclang-5.0-dev; fi
if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'clang++' && $MPI != 'none' ]]; then export CXXFLAGS="-nostdinc++ ${CXXFLAGS}"; fi

# redefine CXX to the actual version used
if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'clang++' ]]; then export CXX=clang++-5.0; fi
if [[ $TRAVIS_OS_NAME == 'linux' && $COMPILER == 'g++'     ]]; then export CXX=g++-6; fi
# downloads and setups local clang on osx
if [[ $TEST_SUITE == 'osx_local_clang' ]]; then . ./.travis_scripts/setup_local_clang.sh; fi

if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then export CXXFLAGS="-DBOOST_HAS_INT128=1 ${CXXFLAGS}"; fi

# cmake 
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then wget https://github.com/Kitware/CMake/releases/download/v3.13.2/cmake-3.13.2-Linux-x86_64.sh; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo sh cmake-3.13.2-Linux-x86_64.sh --prefix=/usr/local --exclude-subdir; fi

# MPI
if [[ $MPI == 'mpich'    ]]; then sudo $apt_get_install mpich libmpich-dev; fi
if [[ $MPI == 'lam'      ]]; then sudo $apt_get_install lam-runtime lam4-dev; fi
if [[ $MPI == 'openmpi'  ]]; then sudo $apt_get_install openmpi-bin libopenmpi-dev; fi
  if [[ $MPI == 'mvapich2' ]]; then 
    ls -A ${DEPS_DIR}/mvapich2-2.3b
    if [[ -z "$(ls -A ${DEPS_DIR}/mvapich2-2.3b)" ]]; then
      wget http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.3b.tar.gz;
      tar xf mvapich2-2.3b.tar.gz;
      cd mvapich2-2.3b;
      if [[ $COMPILER == 'g++' ]]; then ./configure --disable-fortran --enable-cxx --enable-threads=multiple --with-device=ch3:sock CC=gcc-6 CXX=g++-6 --prefix=${DEPS_DIR}/mvapich2-2.3b ; fi 
      if [[ $COMPILER == 'clang++' ]]; then ./configure --disable-fortran --enable-cxx --enable-threads=multiple --with-device=ch3:sock CC=clang-5.0 CXX=clang++-5.0 --prefix=${DEPS_DIR}/mvapich2-2.3b ; fi 
      make -j4;  
      make install;
      cd ..;
    else
      echo "Using cached mvapich2."
    fi
    export PATH=${DEPS_DIR}/mvapich2-2.3b/bin:${PATH}
    # LIBRARY_PATH for clang?osx?
    export LD_LIBRARY_PATH=${DEPS_DIR}/mvapich2-2.3b/lib:${LD_LIBRARY_PATH}
    export LD_RUN_PATH=${DEPS_DIR}/mvapich2-2.3b/lib:${LD_RUN_PATH}
    export LIBRARY_PATH=${DEPS_DIR}/mvapich2-2.3b/lib:${LIBRARY_PATH}
  fi

if [[ $MPI != 'none'    ]]; then export CXX=${DEPS_DIR}/mvapich2-2.3b/bin/mpic++ ; fi # full path, since libtool in hdf5 installation does not understand PATH set above (?)
if [[ $MPI != 'none'    ]]; then export CC=${DEPS_DIR}/mvapich2-2.3b/bin/mpicc ; fi

# no MPI
if [[ $TRAVIS_OS_NAME == 'linux' && $MPI == 'none' ]]; then 
  sudo $apt_get_install boost1.61
  sudo ln -s /usr/lib/x86_64-linux-gnu/libboost_python-py35.so /usr/lib/x86_64-linux-gnu/libboost_python3.so # different naming conventions for boost python with python 3
fi

# for MPI we need boost>=1.59 with mpi support, boost installation based on https://github.com/boostorg/compute/blob/master/.travis.yml
if [[ $TRAVIS_OS_NAME == 'linux' && $MPI != 'none' ]]; then 
  ls -A ${DEPS_DIR}/boost
  if [[ -z "$(ls -A ${DEPS_DIR}/boost)" ]]; then
    wget http://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.tar.gz 
    tar xf boost_1_65_1.tar.gz
    cd boost_1_65_1
    # configure and install
    if [[ $COMPILER == 'g++' ]]; then echo "using gcc : 6.2 : g++-6 ;" > $HOME/user-config.jam; fi
    if [[ $COMPILER == 'clang++' ]]; then echo "using clang : 5.0 : clang++-5.0 ;" > $HOME/user-config.jam; fi
    echo "using mpi : $CC ;" >> $HOME/user-config.jam
    cat $HOME/user-config.jam
    if [[ $COMPILER == 'g++' ]]; then
      ./bootstrap.sh --prefix=${DEPS_DIR}/boost/ --with-libraries=serialization,mpi,thread,date_time,system,iostreams,timer,filesystem
      ./b2 -d0 install
    fi
    if [[ $COMPILER == 'clang++' ]]; then 
      #clang installation taken from https://gist.github.com/jimporter/10442880
      ./bootstrap.sh --prefix=${DEPS_DIR}/boost/ --with-libraries=serialization,mpi,thread,date_time,system,iostreams,timer,filesystem --with-toolset=clang
      ./b2 clean
      ./b2 toolset=clang cxxflags="-std=c++14 -stdlib=libc++" linkflags="-stdlib=libc++" --prefix=${DEPS_DIR}/boost/ -j 4 stage release
      ./b2 install toolset=clang cxxflags="-std=c++14 -stdlib=libc++" linkflags="-stdlib=libc++" --prefix=${DEPS_DIR}/boost/
    fi
    cd ..
  else
    echo "Using cached boost."
  fi
  export BOOST_ROOT=${DEPS_DIR}/boost
  export LD_LIBRARY_PATH=${DEPS_DIR}/boost/lib:${LD_LIBRARY_PATH}
  export LD_RUN_PATH=${DEPS_DIR}/boost/lib:${LD_RUN_PATH}
  export LIBRARY_PATH=${DEPS_DIR}/boost/lib:${LIBRARY_PATH}
  export CPATH=${DEPS_DIR}/boost/include:${CPATH}
fi

# blitz
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew install blitz; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libblitz0-dev; fi

# Ubuntu dependency issue fix
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" libpango-1.0-0 libpangocairo-1.0-0; fi

# hdf5
#if [[ $TRAVIS_OS_NAME == 'linux' && $MPI == 'none' ]]; then sudo $apt_get_install libhdf5-7; fi # Ubuntu dependency issue fix
if [[ $TRAVIS_OS_NAME == 'linux' && $MPI == 'none' ]]; then sudo $apt_get_install libhdf5-serial-dev; fi
if [[ $TRAVIS_OS_NAME == 'linux' && $MPI == 'none' ]]; then sudo $apt_get_install hdf5-tools; fi

# C++ support missing in Debian package ...
#if [[ $TRAVIS_OS_NAME == 'linux' && $MPI != 'none' ]]; then sudo $apt_get_install libhdf5-openmpi-dev; fi 
# ... so we are installing it manually:
  if [[ $TRAVIS_OS_NAME == 'linux' && $MPI != 'none' ]]; then 
    ls -A ${DEPS_DIR}/hdf5
    if [[ -z "$(ls -A ${DEPS_DIR}/hdf5)" ]]; then
      wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar
      tar xf hdf5-1.10.5.tar
      cd hdf5-1.10.5
      CXXFLAGS=-w CFLAGS=-w ./configure --enable-parallel --enable-cxx --enable-unsupported --enable-threadsafe --prefix=${DEPS_DIR}/hdf5/
      make
      sudo make install
      cd ..
    else
      echo "Using cached hdf5."
    fi
    export HDF5_ROOT=${DEPS_DIR}/hdf5
    export LD_LIBRARY_PATH=${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}
    export LD_RUN_PATH=${HDF5_ROOT}/lib:${LD_RUN_PATH}
    export LIBRARY_PATH=${HDF5_ROOT}/lib:${LIBRARY_PATH}
    export CPATH=${HDF5_ROOT}/include:${CPATH}
    export PATH=${HDF5_ROOT}/bin:${PATH}
  fi
