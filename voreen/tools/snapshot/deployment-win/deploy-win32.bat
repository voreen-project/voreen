::------------------------------------------------------
:: Deploys the VoreenVE binary and all dependencies
:: and associated resources to the specified directory.
::------------------------------------------------------

@SET OUT_DIR=%~1
@SET TARGET_MACHINE=%~2
@SET LINKING_MODE=%~3
@SET VS_VERSION=%~4
@SET VS_PATH=%~5
@SET QT_PATH=%~6
@IF "%OUT_DIR%"=="" GOTO usage
@IF "%TARGET_MACHINE%"=="" GOTO usage
@IF NOT "%TARGET_MACHINE%"=="x86" IF NOT "%TARGET_MACHINE%"=="x64" GOTO usage
@IF "%LINKING_MODE%"=="" GOTO usage
@IF NOT "%LINKING_MODE%"=="dynamic-libs" IF NOT "%LINKING_MODE%"=="static-libs" GOTO usage
@IF "%VS_VERSION%"=="" GOTO usage
@IF NOT "%VS_VERSION%"=="VS2008" IF NOT "%VS_VERSION%"=="VS2010" GOTO usage
@IF "%VS_PATH%"=="" GOTO usage
@IF "%QT_PATH%"=="" GOTO usage
@GOTO checkfiles

:usage
@echo.
@echo Error in script usage. The correct usage is:
@echo   %0 OUT_DIR x86^|x64 dynamic-libs^|static-libs VS2008^|VS2010 VS_PATH QT_PATH [module0] [module1] ...
@echo:
@echo For example:
@echo   %0 voreenve-2.6 x86 dynamic-libs VS2008 "C:\Program Files (x86)\Microsoft Visual Studio 9.0" "C:\Qt\2010.04-vc\qt" base plotting zip
@GOTO :EOF


:checkfiles
@IF %VS_VERSION%==VS2008 GOTO checkfiles_VS2008
@IF %VS_VERSION%==VS2010 GOTO checkfiles_VS2010
:checkfiles_VS2008
@IF NOT EXIST "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\Microsoft.VC90.CRT.manifest GOTO vsnotfound
@IF NOT EXIST "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\msvcm90.dll GOTO vsnotfound
@IF NOT EXIST "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\msvcp90.dll GOTO vsnotfound
@IF NOT EXIST "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\msvcr90.dll GOTO vsnotfound
@GOTO checkfiles_Qt
:checkfiles_VS2010
@IF NOT EXIST "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC100.CRT\msvcp100.dll GOTO vsnotfound
@IF NOT EXIST "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC100.CRT\msvcr100.dll GOTO vsnotfound
@GOTO checkfiles_Qt

:checkfiles_Qt
@IF NOT EXIST "%QT_PATH%"\bin\QtCore4.dll GOTO qtnotfound
@IF NOT EXIST "%QT_PATH%"\bin\QtGui4.dll GOTO qtnotfound
@IF NOT EXIST "%QT_PATH%"\bin\QtOpenGL4.dll GOTO qtnotfound
@GOTO deploy

:deploy
@echo Deploying VoreenVE to directory "%OUT_DIR%" ...
@echo 1. Copy deployment files to output directory ...
@SET INCLUDE_FILES=tools/snapshot/deployment-win/deploymentlist.txt
@IF %LINKING_MODE%==dynamic-libs @IF EXIST tools/snapshot/deployment-win/deploymentlist.dynamic.txt @SET INCLUDE_FILES=%INCLUDE_FILES% tools/snapshot/deployment-win/deploymentlist.dynamic.txt 
@IF %LINKING_MODE%==static-libs @IF EXIST tools/snapshot/deployment-win/deploymentlist.static.txt @SET INCLUDE_FILES=%INCLUDE_FILES% tools/snapshot/deployment-win/deploymentlist.static.txt 
@SHIFT
@SHIFT
@SHIFT
:moduleLoop
@IF "%~1"=="" @GOTO loopContinue
   @IF EXIST modules/%~1/deploymentlist.txt @SET INCLUDE_FILES=%INCLUDE_FILES% modules/%~1/deploymentlist.txt
   @IF %TARGET_MACHINE%==x86 @IF EXIST modules/%~1/deploymentlist.win32.txt @SET INCLUDE_FILES=%INCLUDE_FILES% modules/%~1/deploymentlist.win32.txt
   @IF %TARGET_MACHINE%==x64 @IF EXIST modules/%~1/deploymentlist.win64.txt @SET INCLUDE_FILES=%INCLUDE_FILES% modules/%~1/deploymentlist.win64.txt
   @SHIFT
@GOTO moduleLoop
:loopContinue
perl tools/snapshot/copy-files.pl "%OUT_DIR%" %INCLUDE_FILES%
@IF ERRORLEVEL 1 GOTO error

@echo.
@echo 2. Clean output directory "%OUT_DIR%" from excluded files/dirs ...
perl tools/snapshot/clean-directory.pl "%OUT_DIR%" tools/snapshot/snapshot-exclude.txt	
@IF ERRORLEVEL 1 GOTO error

@echo.
@echo 3. Deploying VS runtime library to "%OUT_DIR%" ...
@IF %VS_VERSION%==VS2008 GOTO deployRuntime_VS2008
@IF %VS_VERSION%==VS2010 GOTO deployRuntime_VS2010
:deployRuntime_VS2008
copy "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\Microsoft.VC90.CRT.manifest %OUT_DIR%
copy "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\msvcm90.dll %OUT_DIR%
copy "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\msvcp90.dll %OUT_DIR%
copy "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC90.CRT\msvcr90.dll %OUT_DIR%
@IF ERRORLEVEL 1 GOTO error
@GOTO deployQt
:deployRuntime_VS2010
copy "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC100.CRT\msvcp100.dll %OUT_DIR%
copy "%VS_PATH%"\VC\redist\%TARGET_MACHINE%\Microsoft.VC100.CRT\msvcr100.dll %OUT_DIR%
@IF ERRORLEVEL 1 GOTO error
@GOTO deployQt

:deployQt
@echo.
@echo 4. Deploying Qt libraries to "%OUT_DIR%" ...
copy "%QT_PATH%"\bin\QtCore4.dll %OUT_DIR%
copy "%QT_PATH%"\bin\QtGui4.dll %OUT_DIR%
copy "%QT_PATH%"\bin\QtOpenGL4.dll %OUT_DIR%
@IF ERRORLEVEL 1 GOTO error
mkdir %OUT_DIR%\licenses\Qt
copy "%QT_PATH%"\LICENSE.LGPL %OUT_DIR%\licenses\Qt
::@IF ERRORLEVEL 1 GOTO error

@echo.
@echo Deployed VoreenVE win32-%TARGET_MACHINE% binary to directory "%OUT_DIR%".
@GOTO :EOF

:vsnotfound
@echo.
@echo Visual Studio run time library not found. VS version '%VS_VERSION%' or path '%VS_PATH%' incorrect?
@GOTO :EOF

:qtnotfound
@echo.
@echo Qt libraries not found. Qt path '%QT_PATH%' incorrect?
@GOTO :EOF

:error
@echo.
@echo FAILED to deploy VoreenVE!!!
  