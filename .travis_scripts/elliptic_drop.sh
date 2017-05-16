#!/usr/bin/env sh
set -e
# installing
cd libmpdata++/build
sudo make install
cd ../..
# shallow-water-elliptic-drop
git clone --depth=1 git://github.com/igfuw/shallow-water-elliptic-drop.git
cd shallow-water-elliptic-drop/numerical
mkdir build
cd build
if [[ $TRAVIS_OS_NAME == 'linux' && $CXX == 'clang++' ]]; then 
  cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_BUILD_TYPE=Release ../; 
else
  cmake -DCMAKE_BUILD_TYPE=Release ../; 
fi
make
#- ./spreading_drop_2d_el 
# TODO: make test!
cd ../../..
