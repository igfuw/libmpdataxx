#!/usr/bin/env sh
set -ex
# installing
cd libmpdata++/build
sudo make install
cd ../..
# shallow-water-elliptic-drop
git clone --depth=1 git://github.com/igfuw/shallow-water-elliptic-drop.git
cd shallow-water-elliptic-drop/numerical
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release #TODO: Debug mode?
make
#- ./spreading_drop_2d_el 
# TODO: make test!
cd ../../..
