#!/bin/bash
#cd src || exit

make clean

make

export LD_LIBRARY_PATH=$(pwd)/build:$LD_LIBRARY_PATH

echo 'DONT WORRY!'

./build/fugue

echo 'BE HAPPY!'
