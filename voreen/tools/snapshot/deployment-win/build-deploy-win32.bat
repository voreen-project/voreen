:: ********************************************************************
::  Makes a deployment build of VoreenVE using MSVC 2008/2010 (nmake).
::  Run this script from Voreen root directory. Optionally,
::  this script can used to generate a config.txt 
::  for creating a deployment build with Visual Studio.
::
::  Warning: config.txt is overwritten!
::
::  Parameters:
::    %1: QMAKESPEC: "win32-msvc2008" or "win32-msvc2010"
::    %2: Target machine: "x86" or "x64"
::    %3: Linking mode: "dynamic-libs" or "static-libs"
::    %4: Visual Studio 2008/2010 path
::    %5: Path to qmake.exe
::    %6: "--generate-config" (optional, does only create config.txt)
::    %7-%n: Names of the modules to include in the build. If no module
::           name is passed, the "default" meta-module is used. 
::
:: 
::  Sample call (remove line breaks):
::    call tools\snapshot\deployment-win\build-deploy-win32.bat 
::         win32-msvc2008 x86 dynamic-libs       
::         "C:\Program Files (x86)\Microsoft Visual Studio 9.0"
::         "C:\Qt\2010.04-vc\qt\bin"
::         --generate-config
::         base connectedcomponents devil dicom ffmpeg flowanalysis 
::         fontrendering plotting pvm python segy tiff zip
:: ********************************************************************

:: alias command line params
@SET QMAKESPEC=%~1
@SET TARGET_MACHINE=%~2
@SET LINKING_MODE=%~3
@SET VS_PATH=%~4
@SET QMAKE_PATH=%~5
@SET CONFIG_MODE=0
@IF "%~6"=="--generate-config" @SET CONFIG_MODE=1

:: check params
@IF "%QMAKESPEC%"=="" GOTO usage
@IF NOT "%QMAKESPEC%"=="win32-msvc2008" IF NOT "%QMAKESPEC%"=="win32-msvc2010" GOTO usage
@IF "%TARGET_MACHINE%"=="" GOTO usage
@IF NOT "%TARGET_MACHINE%"=="x86" IF NOT "%TARGET_MACHINE%"=="x64" GOTO usage
@IF "%LINKING_MODE%"=="" GOTO usage
@IF NOT "%LINKING_MODE%"=="dynamic-libs" IF NOT "%LINKING_MODE%"=="static-libs" GOTO usage
@IF "%VS_PATH%"=="" GOTO usage
@IF "%QMAKE_PATH%"=="" GOTO usage
@GOTO checkfiles

:usage
@echo.
@echo Error in script usage. The correct usage is:
@echo   %0 win32-msvc2008^|win32-msvc2010 x86^|x64 dynamic-libs^|static-libs VS_PATH QMAKE_PATH [--generate-config] [module0] [module1] ...
@echo.
@echo For example:
@echo   %0 win32-msvc2008 x86 dynamic-libs "C:\Program Files (x86)\Microsoft Visual Studio 9.0" "C:\Qt\2010.04-vc\qt\bin" default devil dicom
@GOTO :EOF

:checkfiles
@IF NOT EXIST "%VS_PATH%"\VC\vcvarsall.bat GOTO vsnotfound
@IF NOT EXIST "%QMAKE_PATH%"\qmake.exe GOTO qmakenotfound
@GOTO checkconfig

:vsnotfound
@echo.
@echo Visual Studio run time library not found. VS path %VS_PATH% incorrect?
@GOTO :EOF

:qmakenotfound
@echo.
@echo qmake.exe not found. Qmake path %QMAKE_PATH% incorrect?
@GOTO :EOF

:checkconfig
@echo.
@IF NOT EXIST config.txt GOTO :build
@echo WARNING: config.txt will be overwritten. Continue? 
@CHOICE /c:JN
@IF ERRORLEVEL 2 GOTO :EOF

@echo on

:build
:: setup config.txt
copy config-default.txt config.txt
@IF ERRORLEVEL 1 GOTO :EOF

echo VRN_PROJECTS = tgt core qt voreenve >> config.txt
IF %CONFIG_MODE%==0 echo CONFIG += nmake >> config.txt
echo DEFINES += VRN_DEPLOYMENT >> config.txt
echo DEFINES -= VRN_DEBUG >> config.txt
echo DEFINES -= VRN_DYNAMIC_LIBS >> config.txt
echo DEFINES -= VRN_STATIC_LIBS >> config.txt
IF "%LINKING_MODE%" == "dynamic-libs" echo DEFINES += VRN_DYNAMIC_LIBS >> config.txt
IF "%LINKING_MODE%" == "static-libs" echo DEFINES += VRN_STATIC_LIBS >> config.txt

:: modules
@SHIFT 
@SHIFT
@SHIFT
@SHIFT
@SHIFT
@IF "%~1"=="--generate-config" @SHIFT
echo VRN_MODULES = >> config.txt
@IF "%~1"=="" @GOTO defaultModules
:moduleLoop
@IF "%~1"=="" @GOTO continue
   echo VRN_MODULES += "%~1" >> config.txt 
   @SHIFT
@GOTO moduleLoop
:defaultModules
echo VRN_MODULES += default >> config.txt
:continue

:: jump over build instructions, if in config generation mode
@IF %CONFIG_MODE%==1 @GOTO configMode 

:: setup visual studio build environment
call "%VS_PATH%"\VC\vcvarsall.bat %TARGET_MACHINE%
@IF ERRORLEVEL 1 GOTO :EOF

:: run qmake and build with nmake
"%QMAKE_PATH%"\qmake -spec %QMAKESPEC% voreen.pro
@IF ERRORLEVEL 1 GOTO :EOF

nmake clean
nmake release
GOTO :EOF

:configMode
@echo.
@echo -- Generated config.txt for deployment build in Visual Studio. --
