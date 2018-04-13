#!/bin/sh

# Original Repo with manual: https://bitbucket.org/alfonse/glloadgen
# Repo used by voreen with os x patch: https://github.com/mosra/glloadgen

OPENGL_VERSION="4.4"
OPENGL_EXTENSIONS="EXT_texture_compression_s3tc EXT_texture_sRGB EXT_texture_filter_anisotropic"

# Update to newest revision before generating the header
if [ ! -d "glloadgen" ]; then
    echo "Downloading glLoadGen..."
    git clone https://github.com/mosra/glloadgen
else
    echo "Updating glLoadGen..."
    cd glloadgen
    git pull
    cd ..
fi

# generate the header
echo "Generating header..."
lua ./glloadgen/LoadGen.lua -style=pointer_c -spec=gl -version=$OPENGL_VERSION -exts $OPENGL_EXTENSIONS -profile=core core

# rename for c++ compilation
echo "Renaming..."
mv gl_core.c gl_core.cpp

# patch for windows (export functions using TGT_API)
echo "Patching..."
perl -p -i -e 's/(#ifndef APIENTRY)/#include \"tgt\/types.h\"\n\n$1/g' gl_core.h    # add include of tgt/types.h for TGT_API definition
perl -p -i -e 's/(extern.+CODEGEN_FUNCPTR)/TGT_API $1/g' gl_core.h                  # prepend every function with TGT_API

echo "Done."
