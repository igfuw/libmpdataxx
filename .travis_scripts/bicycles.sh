#!/usr/bin/env sh
set -e

# making Python 2 back the default if needed - TODO: support Python3 in libcloudph++
if [[ $PY3DEB != '' ]]; then sudo update-alternatives --remove python /usr/bin/python3; fi
if [[ $PY3DEB != '' ]]; then sudo update-alternatives --install /usr/bin/python python /usr/bin/python2 10; fi

# libcloudph++'s dependencies
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libboost-python-dev python-numpy; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then brew install boost-python; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then git clone --depth=1 git://github.com/thrust/thrust.git; fi
if [[ $TRAVIS_OS_NAME == 'osx' ]]; then sudo ln -s `pwd`/thrust/thrust /usr/local/include/thrust; fi

# odeint
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then git clone --depth=1 https://github.com/boostorg/odeint.git; fi # get boost odeint > 1.58
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo rm -f /usr/include/boost/numeric/odeint.hpp; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo rm -rf /usr/include/boost/numeric/odeint; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo ln -s `pwd`/odeint/include/boost/numeric/odeint.hpp /usr/include/boost/numeric/odeint.hpp; fi
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo ln -s `pwd`/odeint/include/boost/numeric/odeint  /usr/include/boost/numeric/; fi

# newest thrust
git clone --depth=1 git://github.com/thrust/thrust.git;
sudo ln -s `pwd`/thrust/thrust /usr/local/include/thrust;

# libcloudph++ (needed by icicle, skipping tests and CUDA build)
git clone --depth=1 git://github.com/igfuw/libcloudphxx.git
cd libcloudphxx
mkdir build 
cd build
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ../; fi # Travis default is not the packaged     one
if [[ $TRAVIS_OS_NAME == 'osx' && $CXX == 'g++' ]]; then cmake -DCMAKE_CXX_COMPILER=g++-4.8 ../; fi # the one from homebrew
cmake -DCMAKE_BUILD_TYPE=Release ../ 
make 
sudo make install
cd ../..

#bicycles
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then sudo $apt_get_install libboost-program-options-dev; fi
git clone --depth=1 git://github.com/igfuw/bicycles.git
cd bicycles
mkdir build
cd build
cmake ..
make
make test || cat Testing/Temporary/LastTest.log / # "/" intentional! (just to make cat exit with an error code)
cd ../..
set +e # see https://github.com/travis-ci/travis-ci/issues/6522

