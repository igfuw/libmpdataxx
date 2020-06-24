#!/usr/bin/env sh
set -e

if [ $# -ne 1 ]; then
  echo "UWLCM.sh accepts exactly one argument"
  exit 1
fi

# install libmpata++
cd libmpdata++/build
sudo make install
cd ../..

# newest thrust
git clone --depth=1 git://github.com/thrust/thrust.git;
sudo ln -s `pwd`/thrust/thrust /usr/local/include/thrust;

# libcloudph++ 
git clone --depth=1 git://github.com/igfuw/libcloudphxx.git
cd libcloudphxx
mkdir build 
cd build
# RelWithDebInfo = Release with asserts
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../ 
make 
sudo make install
cd ../..

# UWLCM
# if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libboost-program-options1.55-dev; fi
git clone --depth=1 git://github.com/igfuw/UWLCM.git
cd UWLCM
. .travis_scripts/$1.sh
cd ..
set +e # see https://github.com/travis-ci/travis-ci/issues/6522

