#!/bin/bash
# Removes binaries, header files, data files installed by 
# find_orb and optionally the git repositories of its dependencies.

usage () { 
    echo "Usage: /bin/bash UNINSTALL.sh -d {DOWNLOAD_DIR}"
    echo "-d DOWNLOAD_DIR [REQUIRED]"
    echo "    Path to directory where find_orb and its dependencies were downloaded."
    echo "-g REMOVE_GIT [OPTIONAL]"
    echo "    Remove git repositories of find_orb's dependencies."
    exit 1
}

# Add argument catches, if arguments are not recognized run usage function
dir_flag=false
repo_flag=false
while getopts ":d:g" opt; do
  case $opt in
    d) dir="$OPTARG"; dir_flag=true;;
    g) repo_flag=true;;
    \?) echo "Invalid argument -$OPTARG" >&2; usage;
  esac
done

# Make sure the download directory is passed
if [ "$dir_flag" = false ]
then
    usage
fi

BIN=~/bin
INCLUDE=~/include
LIB=~/lib

echo "Uninstalling find_orb."
rm -rf -v ~/.find_orb

# Remove binaries from ~/bin/
array=("astcheck" "find_orb" "fo" "sat_id")
for i in "${array[@]}"
do
    if [ -f $BIN/$i ]; then
        rm -v $BIN/$i
    else
        echo "$i has already been removed."
    fi
done

# Remove header files from ~/include/
array=("afuncs.h" "cgi_func.h" "date.h" "jpleph.h" \
       "mpc_func.h" "showelem.h" "watdefs.h" \
       "brentmin.h" "comets.h" "get_bin.h" \
       "lunar.h" "norad.h" "vislimit.h")
for i in "${array[@]}"
do
    if [ -f $INCLUDE/$i ]; then
        rm -v $INCLUDE/$i
    else
        echo "$i has already been removed."
    fi
done

# Remove supporting libraries
array=("libjpl.a" "liblunar.a" "libsatell.a")
for i in "${array[@]}"
do
    if [ -f $LIB/$i ]; then
        rm -v $LIB/$i
    else
        echo "$i has already been removed."
    fi
done

# If desired, remove git repositories
repos=("lunar" "jpl_eph" "sat_code" "miscell")
if [ "$repo_flag" = true ]; then
    for repo in "${repos[@]}"; do
        if [ -d $dir/$repo ]
        then
            echo "Removing $repo..."
            rm -rf -v $dir/$repo
        else
            echo "$repo repository has already been removed."
        fi
    done
fi
exit