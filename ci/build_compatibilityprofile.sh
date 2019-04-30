#!/bin/bash

# clear build directory (also clears CMake cache)
BUILD_DIR=$1
rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR

# build configuration
build_options=(
    -DCMAKE_BUILD_TYPE=Debug

    -DVRN_OPENGL_COMPATIBILITY_PROFILE=ON #required to build some modules

    -DVRN_BUILD_VOREENVE=ON
    -DVRN_BUILD_VOREENTOOL=ON
    -DVRN_BUILD_SIMPLEGLUT=ON
    -DVRN_BUILD_SIMPLEQT=ON
    -DVRN_BUILD_TESTAPPS=ON
    -DVRN_BUILD_BLASTEST=ON
    -DVRN_BUILD_ITKWRAPPER=ON
    #-DVRN_DEPLOYMENT=ON

    # public modules
    -DVRN_MODULE_BASE=ON
    -DVRN_MODULE_CONNEXE=ON
    -DVRN_MODULE_DYNAMICGLSL=ON #Only compatibility
    -DVRN_MODULE_EXPERIMENTAL=ON #Only compatibility
    -DVRN_MODULE_FLOWREEN=ON
    -DVRN_MODULE_HDF5=ON
    -DVRN_USE_HDF5_VERSION=1.10
    -DVRN_MODULE_PLOTTING=ON
    -DVRN_MODULE_POI=ON
    -DVRN_MODULE_PVM=ON
    -DVRN_MODULE_RANDOMWALKER=ON
    -DVRN_MODULE_SEGY=ON
    -DVRN_MODULE_VOLUMELABELING=ON #Only compatibility
    -DVRN_MODULE_DEVIL=ON
    -DVRN_MODULE_ZIP=ON
    -DVRN_MODULE_FFMPEG=ON
    -DVRN_MODULE_TIFF=ON
    #-DVRN_MODULE_TOUCHTABLE=ON
    -DVRN_MODULE_PFSKEL=OFF
    -DVRN_MODULE_PYTHON=ON
    -DVRN_MODULE_OPENCL=OFF #Until we figure out how to run this in ci
    -DVRN_MODULE_OPENMP=ON
    -DVRN_MODULE_GDCM=ON
    #-DVRN_MODULE_ITK=ON
    -DVRN_MODULE_STAGING=ON
    -DVRN_MODULE_SURFACE=ON
    -DVRN_MODULE_DEPRECATED=ON #Only compatibility
    -DVRN_MODULE_STEREOSCOPY=ON

    # custom modules
    -DVRN_MODULE_GRAPHLAYOUT=ON

    # do not block on failed assertions
    -DVRN_NON_INTERACTIVE=ON
)
cmake "${build_options[@]}" ../voreen

N_CORES=$(nproc)
# start build once again to make errors (if any) more readable
nice make -j$N_CORES || make