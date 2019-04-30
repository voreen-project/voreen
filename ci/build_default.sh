#!/bin/bash

# clear build directory (also clears CMake cache)
BUILD_DIR=$1
rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR

# build configuration
build_options=(
    #just the defaults!
)
cmake "${build_options[@]}" ../voreen

N_CORES=$(nproc)
# start build once again to make errors (if any) more readable
nice make -j$N_CORES || make
