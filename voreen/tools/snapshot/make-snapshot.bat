::------------------------------------------------
:: Creates a source snapshot for win32 or unix
:: in the specified output directory.
:: To be run from voreen root directory.
::
:: NOTE: Make sure to run this script on a clean
::       checkout without any build residues!
::------------------------------------------------

@if "%1" == "" GOTO usage
@IF "%2" == "" GOTO usage
@IF "%1" == "win32" GOTO win32 
@IF "%1" == "unix"  GOTO unix 
@IF NOT "%~3"=="" IF NOT "%~3"=="--force" GOTO usage

:usage
@echo:
@echo Error in script usage. The correct usage is:
@echo     %0 win32^|unix OUT_DIR [--force]
@echo:
@echo For example:
@echo     %0 win32 voreen-src-2.6 
@GOTO :EOF


:win32
@echo.
@echo Creating win32 source snapshot in directory "%2%" ...
@echo.
@echo 1. Copy snapshot files to output directory ...
@IF "%~3"=="--force" perl tools/snapshot/copy-files.pl %2% tools/snapshot/snapshot-include-common.txt tools/snapshot/snapshot-include-win32.txt --force
@IF NOT "%~3"=="--force" perl tools/snapshot/copy-files.pl %2% tools/snapshot/snapshot-include-common.txt tools/snapshot/snapshot-include-win32.txt
@IF ERRORLEVEL 1 GOTO error

@echo 2. Clean output directory "%2%" from excluded files/dirs ...
perl tools/snapshot/clean-directory.pl %2%/include tools/snapshot/snapshot-exclude.txt	
perl tools/snapshot/clean-directory.pl %2%/src tools/snapshot/snapshot-exclude.txt
perl tools/snapshot/clean-directory.pl %2%/apps tools/snapshot/snapshot-exclude.txt
perl tools/snapshot/clean-directory.pl %2%/modules tools/snapshot/snapshot-exclude.txt	
@IF ERRORLEVEL 1 GOTO error

@echo 3. Running convert-eol on snapshot directory ...
perl tools/snapshot/convert-eol.pl %2% win
@IF ERRORLEVEL 1 GOTO error

@echo.
@echo Created win32 source snapshot in directory "%2%".
@GOTO :EOF


:unix
@echo.
@echo Creating unix source snapshot in directory "%2%" ...
@echo.
@echo 1. Copy snapshot files to output directory "%2%" ...
@IF "%~3"=="--force" perl tools/snapshot/copy-files.pl "%2%" tools/snapshot/snapshot-include-common.txt tools/snapshot/snapshot-include-unix.txt --force
@IF NOT "%~3"=="--force" perl tools/snapshot/copy-files.pl "%2%" tools/snapshot/snapshot-include-common.txt tools/snapshot/snapshot-include-unix.txt
@IF ERRORLEVEL 1 GOTO error

@echo 2. Clean output directory from .svn directories ...
perl tools/snapshot/clean-directory.pl %2% tools/snapshot/snapshot-exclude.txt	
@IF ERRORLEVEL 1 GOTO error

@echo 3. Running convert-eol on snapshot directory ...
perl tools/snapshot/convert-eol.pl %2% unix 
@IF ERRORLEVEL 1 GOTO error

@echo.
@echo Created unix source snapshot in directory "%2%".
@GOTO :EOF


:error
@echo.
@echo FAILED to create %1% source snapshot!
  