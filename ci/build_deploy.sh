#!/bin/bash

VRN_DIR=$(pwd)/voreen

INSTALL_DIR=$(realpath $2)
rm -rf $INSTALL_DIR
mkdir $INSTALL_DIR

# clear build directory (also clears CMake cache)
BUILD_DIR=$1
rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR

# build configuration
build_options=(
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
    -DVRN_DEPLOYMENT=ON
    -DVRN_ADD_INSTALL_TARGET=ON

    #set custom module options here
    -DVRN_MODULE_FFMPEG=ON
    -DVRN_MODULE_HDF5=ON
    -DVRN_USE_HDF5_VERSION=1.10
    -DVRN_MODULE_OPENCL=ON
    -DVRN_MODULE_OPENMP=ON
    -DVRN_MODULE_RANDOMWALKER=ON
    -DVRN_MODULE_STEREOSCOPY=ON
    -DVRN_MODULE_TIFF=ON

)
cmake "${build_options[@]}" $VRN_DIR

N_CORES=$(nproc)
# start build once again to make errors (if any) more readable
nice make -j$N_CORES install || make install
