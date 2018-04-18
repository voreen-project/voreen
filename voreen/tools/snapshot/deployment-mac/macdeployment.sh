#!/bin/bash

######################################################################
#                                                                    #
# Voreen - Mac deployment script                                     #
#                                                                    #
# Copyright (C) 2005-2010 Visualization and Computer Graphics Group, #
# Department of Computer Science, University of Muenster, Germany.   #
# <http://viscg.uni-muenster.de>                                     #
#                                                                    #
######################################################################

# Configuration
DEPLOYMENT_TOOL="../../tools/snapshot/deployment-mac/macdeployqt-voreen/macdeployqt-voreen"

APP_BUNDLE_NAME="voreenve.app"
APP_BUNDLE_PATH="."
APP_BINARY_NAME="voreenve"
APP_BUNDLE_DEPLOYMENT_PATH="VoreenVE-2.6.2"

DMG_NAME="VoreenVE-2.6.2"

SHADERS_PATH="../../glsl"



# Delete previous dmg and dist directory
rm $DMG_NAME.dmg 2> /dev/null
rm -rf dist 2> /dev/null



# Deploy app bundle
echo ""
echo "* ---------------------------------------------------"
echo "*  Deploying dependencies of $APP_BUNDLE_NAME ...    "
echo "* ---------------------------------------------------"

# Create temporary distribution directory 'dist' and copy app bundle to it
echo ""
echo "* Creating temporary distribution directory and copying over app bundle ..."
mkdir -v dist
mkdir -v dist/$APP_BUNDLE_DEPLOYMENT_PATH
cp -Rv $APP_BUNDLE_PATH/$APP_BUNDLE_NAME dist/$APP_BUNDLE_DEPLOYMENT_PATH

# Call deployment tool on app bundle (copies dependencies into the app bundle and registers them)
$DEPLOYMENT_TOOL dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME

# deployment tool does currently only consider Qt frameworks and dylibs, but no other frameworks
# => custom frameworks have to be deployed separately
echo ""
echo "* Deploying additional frameworks to dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME ..."
# Python
cp -v /System/Library/Frameworks/Python.framework/Versions/2.5/Python dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Frameworks
install_name_tool -id @executable_path/../Frameworks/python dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Frameworks/python
install_name_tool -change /System/Library/Frameworks/Python.framework/Versions/2.5/Python @executable_path/../Frameworks/python dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/MacOS/$APP_BINARY_NAME


echo ""
echo "* --------------------------------------------------------"
echo "*  Copying over voreen resources to $APP_BUNDLE_NAME ...  "
echo "* --------------------------------------------------------"
#
# Copy readonly resources to app bundle's resources path:
# - shaders
# - textures
# - fonts
# - documentation 
#

mkdir -v dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/textures
rsync -rv --exclude=.svn ../../data/textures/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/textures/
mkdir -v dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/fonts/
rsync -rv --exclude=.svn ../../data/fonts/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/fonts/
mkdir -v dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/doc/
rsync -rv --exclude=.svn ../../doc/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/doc/

#shaders (must be copied to glsl subdir before)
mkdir -v dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/glsl/
rsync -rv --exclude=.svn ../../glsl/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/glsl/

#modules (must be copied to modules subdir before)
mkdir -v dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/modules/
rsync -rv --exclude=.svn ../../modules/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources/modules/

# copy over Info.plist and icon
cp -v ../../tools/macdeployment/Info.plist dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents
cp -v ../../tools/macdeployment/icon.icns dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/Resources

#
# add files from root dir
#
cp ../../Changelog.txt dist/$APP_BUNDLE_DEPLOYMENT_PATH/
cp ../../CREDITS.txt dist/$APP_BUNDLE_DEPLOYMENT_PATH/
cp ../../LICENSE.txt dist/$APP_BUNDLE_DEPLOYMENT_PATH/
cp ../../LICENSE-academic.txt dist/$APP_BUNDLE_DEPLOYMENT_PATH/


#
# add user data to subdir 'data' within deployment path 
#
mkdir -v dist/$APP_BUNDLE_DEPLOYMENT_PATH/data

# copy remaining data
mkdir dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/networks
rsync -rv --exclude=.svn ../../data/networks/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/networks/
mkdir dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/workspaces
rsync -rv --exclude=.svn ../../data/workspaces/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/workspaces/
mkdir dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/scripts
rsync -rv --exclude=.svn ../../data/scripts/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/scripts/
mkdir dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/transferfuncs
rsync -rv --exclude=.svn ../../data/transferfuncs/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/transferfuncs/
mkdir dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/volumes
rsync -rv --exclude=.svn ../../data/volumes/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/volumes/

mkdir dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/misc
rsync -rv --exclude=.svn ../../data/misc/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/data/misc/

#external license files (must be copied to licenses subdir before)
mkdir dist/$APP_BUNDLE_DEPLOYMENT_PATH/licenses
rsync -rv --exclude=.svn ../../licenses/ dist/$APP_BUNDLE_DEPLOYMENT_PATH/licenses/


# create symbolic link to Application directory in disk image
echo ""
echo "* Generating symbolic link to Application directory ..."
ln -sv /Applications dist



# Create disk image
echo ""
echo "* ---------------------------------------------------"
echo "*  Creating $DMG_NAME.dmg from directory 'dist' ...  "
echo "* ---------------------------------------------------"
hdiutil create -volname $DMG_NAME -srcfolder ./dist/ ./$DMG_NAME.dmg


# Finally, print out dylib dependencies
echo ""
echo "* Application dylib dependencies: "
otool -L dist/$APP_BUNDLE_DEPLOYMENT_PATH/$APP_BUNDLE_NAME/Contents/MacOS/$APP_BINARY_NAME


# delete dist directory
#rm -rf dist

echo ""
echo "finished: $DMG_NAME.dmg"
echo ""
