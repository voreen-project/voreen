#!/bin/bash
#################################################
# Creates a source snapshot for win32 or unix.
# The passed snapshot name is used as output 
# directory and as archive name.
#
# To be run from Voreen root directory.
#
# NOTE: Make sure to run this script on a clean
#       checkout without any build residues!
#################################################

function usage() {
  echo -e "\nUsage: $0 <win32|unix|all> <snapshot-name> [--unifdef] [--force]\n"
}

function checkfailure() {
  if [ $1 != 0 ]; then
     echo -e "\nFAILED!\n"
     exit 1
  fi
}

# check input params
if [[ "$1" == "" || "$2" == "" ]]; then
  usage
  exit 1
fi

if [[ $1 != "win32" && $1 != "unix" && $1 != "all" ]]; then
  usage
  exit 1
fi

force=false
forceStr=""
unifdef=false
if [ "$3" == "--force" ]; then
  force=true
  forceStr="--force"
elif [ "$3" == "--unifdef" ]; then
  unifdef=true
elif [ "$3" != "" ]; then
  usage
  exit 1
fi

if [ "$4" == "--force" ]; then
  force=true
  forceStr="--force"
elif [ "$4" == "--unifdef" ]; then
  unifdef=true
elif [ "$4" != "" ]; then
  usage
  exit 1
fi

# check if archive exists
if [ $1 == "win32" ]; then
  archiveName=$2.zip
elif [[ $1 == "unix" || $1 == "all" ]]; then
  archiveName=$2.tar.gz
fi
if [ -e $archiveName ]; then
  if [ $force == false ]; then
    echo -e "\nArchive '$archiveName' already exists. Overwrite? (y=yes, c=cancel) ";
    read userinput;
    if [[ $userinput != "y" && $userinput != "yes" ]]; then 
      echo -e "Cancel.\n"
      exit 1
    fi
  fi
  rm $archiveName
fi


# 1. run copy script
if [ $1 == "win32" ]; then
  echo -e "\nCreating win32 source snapshot in directory '$2' ...\n"
  echo -e "1. Copy snapshot files to output directory ..."
  perl tools/snapshot/copy-files.pl $2 tools/snapshot/snapshot-include-common.txt tools/snapshot/snapshot-include-win32.txt $forceStr
elif [ $1 == "unix" ]; then
  echo -e "\nCreating unix source snapshot in directory '$2' ...\n"
  echo -e "1. Copy snapshot files to output directory ..."
  perl tools/snapshot/copy-files.pl $2 tools/snapshot/snapshot-include-common.txt tools/snapshot/snapshot-include-unix.txt $forceStr
elif [ $1 == "all" ]; then
  echo -e "\nCreating full source snapshot (option 'all') in directory '$2' ...\n"
  echo -e "1. Copy snapshot files to output directory ..."
  perl tools/snapshot/copy-files.pl $2 tools/snapshot/snapshot-include-common.txt tools/snapshot/snapshot-include-unix.txt tools/snapshot/snapshot-include-win32.txt $forceStr
fi
checkfailure $?

# 2. clean output directory
echo -e "\n2. Clean snapshot directory from excluded files and directories ..."
perl tools/snapshot/clean-directory.pl $2 tools/snapshot/snapshot-exclude.txt
checkfailure $?

# 3. run unifdef on output directory
if [ $unifdef == true ]; then
    echo -e "\n3. Running unifdef on snapshot directory ..."
    perl tools/snapshot/unifdef.pl $2/include tools/snapshot/snapshot-unifdef-symbols.txt 
    perl tools/snapshot/unifdef.pl $2/src tools/snapshot/snapshot-unifdef-symbols.txt
    perl tools/snapshot/unifdef.pl $2/apps tools/snapshot/snapshot-unifdef-symbols.txt
    perl tools/snapshot/unifdef.pl $2/modules tools/snapshot/snapshot-unifdef-symbols.txt
    checkfailure $?
else
    echo -e "\n3. Skipping unifdef."
fi

# 4. run convert-eol on output directory
echo -e "\n4. Running convert-eol on snapshot directory ..."
if [ $1 == "win32" ]; then
    perl tools/snapshot/convert-eol.pl $2 "win" 
elif [[ $1 == "unix" || $1 == "all" ]]; then
    perl tools/snapshot/convert-eol.pl $2 "unix" 
fi
checkfailure $?

# 5. create snapshot archive from output directory 
echo -e "\n5. Create snapshot archive '$archiveName' from directory '$2' ..."
if [ $1 == "win32" ]; then
  zip $archiveName -r $2
elif [[ $1 == "unix" || $1 == "all" ]]; then
  tar -czvf $archiveName $2
fi
checkfailure $?

echo -e "\nCreated $1 source snapshot in directory '$2' and archive '$archiveName'.\n"
exit 0

