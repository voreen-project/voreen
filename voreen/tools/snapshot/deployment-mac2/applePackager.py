import subprocess
import sys
import os
import shutil
import re

# Configuration

# only dependencies matching these string are copied to the app bundle (this is a case sensitive substring matching)
#searchDeps = set(["libboost", "Qt", "libtgt", "libvoreen_core", "libvoreen_qt", "libIL", "libicu"])
searchDeps = set(["libtgt", "libvoreen_core", "libvoreen_qt", "/usr/local/"])

########################################################################################################
        
def getQtPlugins(pluginRoot) :
    for root, subdirs, files in os.walk(pluginRoot) :
        for file in files :
            if file.split(".")[-1] == "dylib" :
                yield os.path.join(root,file)
        

def getDependencies(file) :
    deps = []
    o = subprocess.Popen(["/usr/bin/xcrun", "otool", "-L", file], stdout=subprocess.PIPE)
    for l in o.stdout:
        if l[0] == '\t':
            deps.append(l.split(' ', 1)[0][1:])
    return deps

def changeDependenyPath(file, oldPath, newPath) :
    return subprocess.call(["/usr/bin/xcrun", "install_name_tool", "-change", oldPath, newPath, file])
    
    
def copyDependencies(binary, libdir) :
    return copyDependenciesHelper(binary, libdir, set())
    
def copyDependenciesHelper(binary, libdir, deplist) :
    binpath = re.match("(.+/)", binary).group(1)
    for dep in getDependencies(binary) :
        dep = re.sub("@loader_path/", binpath, dep)
        for searchDep in searchDeps :
            if searchDep in dep :
                newDepPath = libDir + dep.split("/")[-1]
                if (not newDepPath in deplist) :
                    shutil.copy(dep, libDir)
                    deplist.add(newDepPath)
                    copyDependenciesHelper(dep, libdir, deplist)
    return deplist
    
def patchDependencies(file) :
    for dep in getDependencies(file) :
            for searchDep in searchDeps :
                if searchDep in dep :
                    changeDependenyPath(file, dep, "@executable_path/../Libraries/" + dep.split("/")[-1])

def patchId(file) :
    # make file writable for patching
    subprocess.call(["chmod", "644", file]);
    subprocess.call(["/usr/bin/xcrun", "install_name_tool", "-id", file.split("/")[-1], file])

def checkDeploymentMode(applicationSourcePath) :
    applicationSource = open(applicationSourcePath)
    deploymentOn = False
    for line in applicationSource :
        if re.match("#define VRN_DEPLOYMENT", line) :
            deploymentOn = True
            break
    return deploymentOn

def enableDeploymentMode(applicationSourcePath) :
    os.system("perl -p -i -e 's/\/\/(#define VRN_DEPLOYMENT$)/$1/m' " + applicationSourcePath)

def getActivatedModules() :
    moduleRegistration = xcodedir + "/gen_moduleregistration.h"
    modules = []
    reg = open(moduleRegistration)
    for line in reg :
        match = re.search("new .+?Module\(\"(.+?)\"\)", line)
        if match :
            modules.append(match.group(1))
    return modules


if len(sys.argv) != 3 :
    print ("usage:\n\tpython applePackager.py xcodeDir qt_root\n")
    print("\txcodeDir is the directory containing xcode project")
    print("\tqt_root is the qt install directory (for example /usr/local/opt/qt)")
    print("")
    exit(1)


# the base path where are the voreen binaries to be packaged
basePath = "../../../"
binaryBasePath = basePath + "bin/Release/"
voreenExecutable = "voreenve"
executablePath = binaryBasePath + voreenExecutable

# the dir where the qt installation is located
xcodedir = sys.argv[1]
qtroot = sys.argv[2]
          
if not qtroot.endswith("/") :
    qtroot = qtroot + "/"
qtpluginDir = qtroot + "plugins"

# check for existence of xcode file
xcodeproj = xcodedir + "/Voreen.xcodeproj"
if not os.path.isdir(xcodeproj) :
    print("Can't find project file: " + xcodeproj)
    exit(1)

# check for existence of qt plugin dir
if not os.path.isdir(qtpluginDir) :
    print "Can't find Qt plugin directory: " + qtpluginDir
    exit(1)

# build the project with deployment setting and release mode
applicationSourcePath = basePath + "src/core/voreenapplication.cpp"
if not checkDeploymentMode(applicationSourcePath) :
    enableDeploymentMode(applicationSourcePath)

if checkDeploymentMode(applicationSourcePath) :
        command = "/bin/sh -c 'cd " + xcodedir + "; /usr/bin/xcrun xcodebuild -scheme voreenve -configuration Release'"
        os.system(command)
else :
        print("VRN_DEPLOYMENT in voreenapplication.cpp could not be activated")
        exit(1)


# check for existence of the executable
if not os.path.isfile(executablePath) :
    print "Can't find %s" % (executablePath)
    exit(1)
    
# define output pathes inside the bundle
outputName = "VoreenVE"
outputPath =  outputName + ".app/"
executableDir = outputPath + "Contents/MacOS/";
libDir = outputPath + "Contents/Libraries/";
resourceDir =  outputPath + "Contents/Resources/";
    
# remove old app bundle and copy skel
skelBundle = "VoreenVE.skel"
shutil.rmtree(outputPath, ignore_errors=1)
shutil.copytree(skelBundle, outputPath)
os.mkdir(executableDir)
os.mkdir(libDir)

# copy voreenve executable  to app bundle
shutil.copy(executablePath, executableDir)
 
# copy all dependencies of executable
deps = copyDependencies(executablePath, libDir)

# patch executable dependencies
patchDependencies(executableDir + voreenExecutable)
    
# patch pathes for dylibs in all voreen binaries
for file in deps :
    patchId(file)
    patchDependencies(file)
    
# copy and patch QtPlugins
shutil.copytree(qtpluginDir, outputPath + "Contents/QtPlugins")
for plugin in getQtPlugins(outputPath + "Contents/QtPlugins/") :
    patchId(plugin)
    patchDependencies(plugin)
    
# copy voreen root
voreenRoot = resourceDir + "voreenRoot/"
shutil.copytree(basePath + "ext/tgt", voreenRoot + "ext/tgt")
shutil.copytree(basePath + "resource", voreenRoot + "resource")
shutil.copytree(basePath + "src", voreenRoot + "src")

# copy only activated modules
modules = getActivatedModules()
for module in modules :
    shutil.copytree(basePath + module, voreenRoot + module)

# create dmg file
os.system("hdiutil create %s.dmg -srcfolder %s -format UDBZ -ov" % (outputName, outputPath))

# remove app bundle
shutil.rmtree(outputPath, ignore_errors=1)

