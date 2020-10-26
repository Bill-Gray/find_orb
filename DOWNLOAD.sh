#!/bin/bash
# This script will download find_orb's dependencies into DOWNLOAD_DIR.
# DOWNLOAD_DIR should be the same directory where the find_orb
# repository was downloaded such that:
# ls DOWNLOAD_DIR
#   find_orb
#
# After running this script DOWNLOAD_DIR should look as follows:
# ls DOWNLOAD_DIR
#   find_orb
#   jpl_eph
#   lunar
#   miscell
#   sat_code

usage () { 
    echo "Usage: /bin/bash DOWNLOAD.sh -d {DOWNLOAD_DIR}"
    echo "-d DOWNLOAD_DIR [REQUIRED]"
    echo "    Path to directory in which to download find_orb's dependencies."
    echo "    Should be the same top-level directory where find_orb repository was downloaded."
    exit 1
}

# Add argument catches, if arguments are not recognized run usage function
dir_flag=false
while getopts ":d:" opt; do
  case $opt in
    d) dir="$OPTARG"; dir_flag=true;;
    \?) echo "Invalid argument -$OPTARG" >&2; usage;
  esac
done

# Make sure the download directory is passed
if [ "$dir_flag" = false ]
then
    usage
fi

# Create list of repositories
repos=("lunar" "jpl_eph" "sat_code" "miscell")

# If the repositories have not been downloaded, download them. 
echo "Cloning required repositories."
for repo in "${repos[@]}"; do
    if [ ! -d ../$repo ]
    then
        git clone https://github.com/Bill-Gray/$repo.git $dir/$repo 
    else
        echo "$repo repository already exists in $dir"
    fi
done
exit