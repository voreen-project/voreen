#!/bin/bash

# Script to generate an ASCII logo using figlet as a c/c++ compatible multi-line
# const char*. The resulting file can be included in the c/c++ source.
#
# Usage:
#
#    ./genlogo.sh figlet_font string outfile
#
# Use quotes to pass multiple words as target string, such as "voreen 4.4".

if [ $# -ne 3 ]
then
    echo "Usage:   ./`basename $0` figlet_font string outfile"
    echo "Example: ./`basename $0` standard voreen logo.h"
    exit 1
fi

BINCALL="figlet -k -f"
CMD="$BINCALL $1 $2 | sed 's/\\/\\\\/g' | sed 's/.$/&\\n\"/' | sed 's/^./\"&/' > $3"
CLEANCMD="$BINCALL $1 $2 >> $3"

# do cleanup of command line so it can be evaluated since
# bash strips single backslashes during eval and echo
BASHCMD=$(echo -E "$CMD" | sed 's/\\/\\\\/g' | sed 's/\"/\\\"/g')
CLEANBASHCMD=$(echo -E "$CLEANCMD" | sed 's/\\/\\\\/g' | sed 's/\"/\\\"/g')

echo "Using $BASHCMD ..."
eval $BASHCMD
echo "//Generated using \"$BASHCMD\". Actual appearance during output:" >> $3
echo "/*" >> $3
eval $CLEANBASHCMD
echo "*/" >> $3
