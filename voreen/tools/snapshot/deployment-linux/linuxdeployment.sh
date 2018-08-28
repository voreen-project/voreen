#!/bin/sh

######################################################################
#                                                                    #
# Voreen - Linux deployment script                                   #
#                                                                    #
# Copyright (C) 2005-2018 University of Muenster, Germany,           #
# Department of Computer Science.                                    #
#                                                                    #
######################################################################

# Configuration
INSTALL_DIR=$1
TOOLS_DIR=$2
DEPLOYMENT_TOOL="linuxdeployqt-continuous-x86_64.AppImage"
DEPLOY_DIR="VoreenVE"

# Prepare deployment
echo ""
echo "* ---------------------------------------------------"
echo "*  Preparing deployment ...                          "
echo "* ---------------------------------------------------"

# Create file tree
cd $TOOLS_DIR
mkdir -v $DEPLOY_DIR
mkdir -vp tmp/usr/bin
mkdir -vp tmp/usr/lib
mkdir -vp tmp/usr/share/applications
mkdir -vp tmp/usr/share/icons/hicolor/256x256

# Copy needed files.
find $INSTALL_DIR -name '*.so' -exec cp {} tmp/usr/lib \;
cp "$INSTALL_DIR/voreenve" tmp/usr/bin/voreenve
cp voreenve.desktop tmp/usr/share/applications/voreenve.desktop
cp voreenve.png tmp/usr/share/icons/hicolor/256x256/voreenve.png

# Prepare deployment
echo ""
echo "* ---------------------------------------------------"
echo "*  Executing deployment tool ...                     "
echo "* ---------------------------------------------------"
# The arch sometimes fails to detect. We are not planning to distribute binaries for other arches anyway
ARCH="x86_64" ./$DEPLOYMENT_TOOL ./tmp/usr/share/applications/voreenve.desktop -appimage -bundle-non-qt-libs -extra-plugins=iconengines/libqsvgicon.so

# Copy resources.
echo ""
echo "* ---------------------------------------------------"
echo "*  Copying resources ...                             "
echo "* ---------------------------------------------------"
mv VoreenVE-x86_64.AppImage $DEPLOY_DIR/VoreenVE-x86_64.AppImage
find $INSTALL_DIR -mindepth 1 -maxdepth 1 -type d | xargs cp -rt $DEPLOY_DIR
find $INSTALL_DIR -maxdepth 1 -type f -exec grep -Iq . {} \; -and -print | xargs cp -t $DEPLOY_DIR

# Delete temp files
rm -rf tmp

# Create archieve
tar -czvf $INSTALL_DIR/$DEPLOY_DIR.tar.gz $DEPLOY_DIR

# Delete
rm -rf $DEPLOY_DIR

echo ""
echo "* ---------------------------------------------------"
echo "*  Done.                                             "
echo "* ---------------------------------------------------"
