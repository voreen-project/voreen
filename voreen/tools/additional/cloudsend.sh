#!/usr/bin/env bash

############################################################
#
# cloudsend.sh
#
# Uses curl to send files to a shared
# Nextcloud/Owncloud folder
#
# Usage: ./cloudsend.sh <file> <folderLink>
# Help:  ./cloudsend.sh -h
#
# Gustavo Arnosti Neves
# https://github.com/tavinus
#
# Get this script to current folder with:
# curl -O 'https://gist.githubusercontent.com/tavinus/93bdbc051728748787dc22a58dfe58d8/raw/cloudsend.sh' && chmod +x cloudsend.sh
#
############################################################

# Modified from original by Dominik Drees


CS_VERSION="0.0.7"

CLOUDURL=""
FOLDERTOKEN=""

PUBSUFFIX="public.php/webdav"
HEADER='X-Requested-With: XMLHttpRequest'
INSECURE=''

# https://cloud.mydomain.net/s/fLDzToZF4MLvG28
# curl -k -T myFile.ext -u "fLDzToZF4MLvG28:" -H 'X-Requested-With: XMLHttpRequest' https://cloud.mydomain.net/public.php/webdav/myFile.ext

printVersion() {
        printf "%s\n" "CloudSender v$CS_VERSION"
}

initError() {
        printVersion >&2
        printf "%s\n" "Init Error! $1" >&2
        printf "%s\n" "Try: $0 --help" >&2
        exit 1
}

usage() {
        printVersion
        printf "\n%s%s\n" "Parameters:" "
  -h | --help      Print this help and exits
  -V | --version   Prints version and exits
  -k | --insecure  Uses curl with -k option (https insecure)"
        printf "\n%s\n%s\n" "Use:" "  $0 <filepath> <folderLink> <password>"
        printf "\n%s\n%s\n\n" "Example:" "  $0 './myfile.txt' 'https://cloud.mydomain.net/s/fLDzToZF4MLvG28' foobar123"
        exit 0
}

##########################
# Process parameters

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        usage
fi

if [ "$1" = "-V" ] || [ "$1" = "--version" ]; then
        printVersion
        exit 0
fi

if [ "$1" = "-k" ] || [ "$1" = "--insecure" ]; then
        INSECURE=' -k'
        printf "%s\n" " > Insecure mode ON"
        shift
fi


##########################
# Validate input

FILENAME="$1"
CLOUDURL="${2%/s/*}"
FOLDERTOKEN="${2##*/s/}"
PASSWORD="$3"
MORE_OPTIONS=""

if [ -d "$FILENAME" ]; then
        MORE_OPTIONS="$MORE_OPTIONS -X MKCOL"
fi

if [ -z "$CLOUDURL" ]; then
        initError "Empty URL! Nowhere to send..."
fi

if [ -z "$FOLDERTOKEN" ]; then
        initError "Empty Folder Token! Nowhere to send..."
fi


##########################
# Check for curl

CURLBIN='/usr/bin/curl'
if [ ! -x "$CURLBIN" ]; then
        CURLBIN="$(which curl 2>/dev/null)"
        if [ ! -x "$CURLBIN" ]; then
                initError "No curl found on system!"
        fi
fi


##########################
# Send file

"$CURLBIN"$INSECURE -T "$FILENAME" -u "$FOLDERTOKEN":"$PASSWORD" -H "$HEADER" "$CLOUDURL/$PUBSUFFIX/$FILENAME" $MORE_OPTIONS
