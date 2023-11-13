#!/usr/bin/env bash
################################################################################
##
## Create release for BISCUIT
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

# Colors for text
# Taken from: https://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
DARK_GRAY='\033[1;30m'
RED='\033[0;31m'
LIGHT_RED='\033[1;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NO_COLOR='\033[0m'

# Where biscuit version is defined
VERSION_FILE="src/biscuit.h"

# Script name
SCRIPT_NAME="$0"

# Helper functions
#-----------------------------------------------------------------------------------------------------------------------
# Print helpful usage information
function usage {
    >&2 echo -e
    >&2 echo -e "Usage: $0 [-h,--help] TYPE"
    >&2 echo -e
    >&2 echo -e "Required inputs:"
    >&2 echo -e "    TYPE [patch|minor|major] : which version number to increment"
    >&2 echo -e
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "    -h,--help : Print help message and exit"
    >&2 echo -e
    >&2 echo -e "Please run from the top directory of project"
    >&2 echo -e
}

# Print green text
function make_green {
    echo -e "${GREEN}$@${NO_COLOR}"
}

# Print light red text
function make_light_red {
    echo -e "${LIGHT_RED}$@${NO_COLOR}"
}

# Print red text
function make_red {
    echo -e "${RED}$@${NO_COLOR}"
}

# Print dark gray text
function make_dark_gray {
    echo -e "${DARK_GRAY}$@${NO_COLOR}"
}

# Print out log message
function message {
    echo -e $(make_dark_gray "[${SCRIPT_NAME}::M]") $@
}

# Print out warning message
function warning {
    echo -e $(make_dark_gray "[${SCRIPT_NAME}::W]") $@
}

# Kill with error message
function die {
    >&2 echo -e $(make_dark_gray "[${SCRIPT_NAME}::E]") $(make_red "$@")
    exit 1
}

# Check version definition file exists
function check_file {
    [ -f ${VERSION_FILE} ] || die "${VERSION_FILE} could not be found!"
}

# Check command line type input
function check_type {
    if [ "$#" -ne "1" ]; then
        usage
        die "Too many arguments given"
    fi

    case "$1" in
        patch )
            ;;
        minor )
            ;;
        major )
            ;;
        * )
            die "Unknown type: $1"
            ;;
    esac
}

function get_version_piece {
    set -- "$1" "$(grep -m 1 "^#define\([ \t]*\)$1" ${VERSION_FILE} | cut -d" " -f3 | sed 's/\"//g')"

    if [[ $1 != "BISCUIT_VERSION_DEVEL" ]] && [ -z "$2" ]; then
        die "Could not retrieve $1 from ${VERSION_FILE}"
    fi

    echo "$2"
}

function get_version {
    VERSION_MAJOR=$(get_version_piece "BISCUIT_VERSION_MAJOR")
    VERSION_MINOR=$(get_version_piece "BISCUIT_VERSION_MINOR")
    VERSION_PATCH=$(get_version_piece "BISCUIT_VERSION_PATCH")
    VERSION_DEVEL=$(get_version_piece "BISCUIT_VERSION_DEVEL")
}

function print_version {
    echo -e "version ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}${VERSION_DEVEL}"
}

function bump_version {
    case "$1" in
        major )
            VERSION_MAJOR=$((VERSION_MAJOR+1))
            VERSION_MINOR=0
            VERSION_PATCH=0
            VERSION_DEVEL=""
            ;;
        minor )
            VERSION_MINOR=$((VERSION_MINOR+1))
            VERSION_PATCH=0
            VERSION_DEVEL=""
            ;;
        patch )
            if [[ ${VERSION_DEVEL} == "" ]]; then
                VERSION_PATCH=$((VERSION_PATCH+1))
            else
                VERSION_DEVEL=""
            fi
            ;;
        dev )
            VERSION_PATCH=$((VERSION_PATCH+1))
            VERSION_DEVEL="-dev"
            ;;
        * )
            die "Unknown type to bump version ($1). Must be either major, minor, patch, or dev."
            ;;
    esac
}

function update_version_piece {
    # Check for inputs
    if [ -z $1 ]; then
        return 0
    fi

    # Set up OLD entry to update to NEW
    CUR=$(get_version_piece $1)
    STR=$2

    # DEVEL entry needs to have some double quotes escaped
    if [[ $1 == "BISCUIT_VERSION_DEVEL" ]]; then
        CUR="\"$(get_version_piece $1)\""
        if [[ -z ${STR} ]]; then
            STR="\"\""
        else
            STR="\"-dev\""
        fi
    fi

    OLD="#define $1 ${CUR}"
    NEW="#define $1 ${STR}"

    # Update version in place
    sed -i "s/${OLD}/${NEW}/g" ${VERSION_FILE} || die "Could not update $1 in ${VERSION_FILE}"

    # Check new version
    if [ -z $2 ]; then
        [ -z $(get_version_piece $1) ] || die "Update version failed for $1"
    else
        [[ $(get_version_piece $1) == "$2" ]] || die "Update version failed for $1"
    fi

    return 0
}

function update_version_file {
    update_version_piece "BISCUIT_VERSION_MAJOR" $1
    update_version_piece "BISCUIT_VERSION_MINOR" $2
    update_version_piece "BISCUIT_VERSION_PATCH" $3
    update_version_piece "BISCUIT_VERSION_DEVEL" $4
}
#-----------------------------------------------------------------------------------------------------------------------

# Main
#-----------------------------------------------------------------------------------------------------------------------
# If no inputs, print usage
if [ "$#" -eq "0" ]; then
    usage
    exit 0
fi

# Verify location script is run from and get current version
check_file
get_version

# Get and check TYPE input
TYPE=""
while [ "$#" -ne "0" ]; do
    case "$1" in
        -h|--help )
            usage
            exit 0
            ;;
        -- )
            shift
            break
            ;;
        -*|--*= )
            die "Unknown input: $1"
            ;;
        * )
            TYPE="${TYPE} $1"
            shift;;
    esac
done

check_type ${TYPE}

# Print current version and new version number
message $(make_green "Current version:") $(print_version)
bump_version ${TYPE}
message $(make_green "Release type:") ${TYPE}
message $(make_green "New version to create:") $(print_version)

# Make sure user approves of the new release number
echo
read -p $'\e[1;33mProceed with creating '"`print_version`"$'? [n/Y]:\e[0m ' go_ahead
if [[ ${go_ahead} != "Y" ]]; then
    warning $(make_light_red "Stopping release creation")
    exit 0
fi

# Set new version in VERSION_FILE
update_version_file ${VERSION_MAJOR} ${VERSION_MINOR} ${VERSION_PATCH} ${VERSION_DEVEL}

# Commit new version
MSG=`print_version`
message $(make_green "Creating commit with message:") ${MSG}
git add ${VERSION_FILE}
git commit -m "${MSG}"

# Create tag for new version (x.y.z.YYYYMMDD)
TAG="v${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}.`date +%Y%m%d`"
message $(make_green "Creating tag:") ${TAG}
git tag "${TAG}"

# Create release zip file with extra git files/dirs removed
RELEASE_DIR=biscuit-release
RELEASE_ZIP=release-source.zip
message $(make_green "Creating zip file")
rm -rf ${RELEASE_DIR} ${RELEASE_ZIP}
git clone . ${RELEASE_DIR}
rm -rf ${RELEASE_DIR}/.git ${RELEASE_DIR}/.github ${RELEASE_DIR}/.gitignore
zip -r ${RELEASE_ZIP} ${RELEASE_DIR}

# Create linux executable
message $(make_green "Creating linux executable")
cd ${RELEASE_DIR}
mkdir build && cd build && cmake ../ && make
cd ../../
rsync ${RELEASE_DIR}/build/src/biscuit biscuit_${VERSION_MAJOR}_${VERSION_MINOR}_${VERSION_PATCH}_linux_amd64

# Update version to devel
bump_version "dev"
message $(make_green "New dev version to create:") $(print_version)

# Set new dev version in VERSION_FILE
update_version_file ${VERSION_MAJOR} ${VERSION_MINOR} ${VERSION_PATCH} ${VERSION_DEVEL}

# Commit devel version
MSG=`print_version`
message $(make_green "Creating commit with message:") ${MSG}
git add ${VERSION_FILE}
git commit -m "${MSG}"
#-----------------------------------------------------------------------------------------------------------------------
