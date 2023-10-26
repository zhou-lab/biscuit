#!/usr/bin/env bash
################################################################################
##
## Check file matches specified sha256 hash
##
## Created by:
##   Jacob Morrison
##
## Creation date:
##   October 2023
##
## Update notes:
##   Oct 2023 -
##     - Initial creation
##
################################################################################

# Check command exists
function this_exists() {
    which "${1}" > /dev/null 2>&1
}

# Print usage information
function usage {
    >&2 echo -e "\nUsage: confirm_download.sh [-h,--help] <filename> <sha_256_hash>\n"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\tfilename     : name of file to check hash"
    >&2 echo -e "\tsha_256_hash : hash to compare against one found for filename\n"
}

# Process command line arguments
OPTS=$(getopt \
    --options h \
    --long help \
    --name "$(basename "$0")" \
    -- "$@"
)
eval set -- ${OPTS}

while true; do
    case "$1" in 
        -h|--help )
            usage
            exit 0
            ;;
        -- )
            shift
            break
            ;;
        * )
            >&2 echo "Unknown option : $1"
            usage
            exit 1
            ;;
    esac
done

# Make sure there are the correct number of inputs
if [[ $# -ne 2 ]]; then
    >&2 echo "$0: Missing inputs"
    usage
    exit 1
fi

# Fill required positional arguments
filename=${1}
shahash=${2}

# Find which sha hash command we can use
sha_command=""
if this_exists sha256sum; then
    sha_command="sha256sum"
elif this_exists shasum; then
    sha_command="shasum -a256"
else
    unset sha_command
fi

if [ -z "${sha_command-}" ]; then
    echo "Could not find shasum command. Skipping verification for ${filename}."
else
    echo "${shahash} ${filename}" | ${sha_command} -c - || { echo "${filename} did not match expected hash. Exiting."; exit 1; }
fi
