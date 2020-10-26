#!/bin/bash
# This script will build and compile find_orb and its dependencies.

usage () { 
    echo "Usage: /bin/bash INSTALL.sh -d {DOWNLOAD_DIR}"
    echo "-d DOWNLOAD_DIR [REQUIRED]"
    echo "    Path to directory where find_orb and its dependencies were downloaded."
    echo "-u UPDATE [OPTIONAL]"
    echo "    Check find_orb's dependencies for updates and pull them."
    exit 1
}

# Add argument catches, if arguments are not recognized run usage function
dir_flag=false
update_flag=false
while getopts ":d:u" opt; do
  case $opt in
    d) dir="$OPTARG"; dir_flag=true;;
    u) update_flag=true;;
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

# Make sure the supporting repos have been downloaded.
# If they have and the update flag has been triggered, 
# then pull the latest updates
for repo in "${repos[@]}"; do
    if [ ! -d $dir/$repo ]
    then
        echo "Could not find the $repo repository!"
        exit 1
    else
        if [ "$update_flag" = true ]
        then
            cd $dir/$repo && git pull origin master
        fi
    fi
done

echo "Building find_orb and its dependencies."
# Make each dependency and find_orb

cd $dir/lunar \
    && make clean \
    && make \
    && make install 
cd $dir/jpl_eph \
    && make clean \
    && make libjpl.a \
    && make install 
cd $dir/lunar \
    && make integrat 
cd $dir/sat_code \
    && make clean \
    && make sat_id \
    && make install 
cd $dir/find_orb \
    && make clean \
    && make \
    && make install 

# Download DE430 if it hasn't already been downloaded.
JPL_EPH_FILE=~/.find_orb/linux_p1550p2650.430t
if [ ! -f "$JPL_EPH_FILE" ]; then 
    cd ~/.find_orb \
        && wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430t/linux_p1550p2650.430t
fi

# Add executables to PATH
export PATH="$PATH:~/bin"
exit 