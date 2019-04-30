#!/bin/bash

# clear build directory (also clears CMake cache)
BUILD_DIR=$1
rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR

# build configuration
build_options=(
    -DVRN_BUILD_VOREENVE=ON
    -DVRN_BUILD_VOREENTOOL=ON
    -DVRN_BUILD_SIMPLEGLUT=ON
    -DVRN_BUILD_SIMPLEQT=ON
    -DVRN_BUILD_TESTAPPS=ON

    -DVRN_SHARED_LIBS=OFF
)
cmake "${build_options[@]}" ../voreen

N_CORES=$(nproc)
# start build once again to make errors (if any) more readable
nice make -j$N_CORES || make
