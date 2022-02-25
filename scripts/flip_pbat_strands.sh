#!/bin/bash
################################################################################
##
## PBAT-style WGBS protocols tend to have reads 1 and 2 flipped relative to
## other protocols. At times, it may be necessary to flip which strand flag the
## reads in the aligned BAM were given. This script is an easy way to do that
## using GNU awk.
##
## Notes:
##   1.) samtools and awk must be in PATH for script to work
##
## Created by:
##   Jacob Morrison
##
## Creation date:
##   January 2021
##
## Update notes:
##
################################################################################

set -euo pipefail

# Check for samtools, awk in PATH
function check_path {
  if [[ `which samtools 2>&1 > /dev/null` ]]; then
      >&2 echo "samtools does not exist in PATH"
      exit 1
  else
      >&2 echo "Using samtools found at: `which samtools`"
  fi
  if [[ `which awk 2>&1 > /dev/null` ]]; then
      >&2 echo "awk does not exist in PATH"
      exit 1
  else
      if awk --version | grep -q GNU; then
          >&2 echo "Using GNU awk found at: `which awk`"
      else
          >&2 echo "It doesn't appear you are using GNU awk"
          >&2 echo "Try adding GNU awk at the front of PATH"
          exit 1
      fi
  fi
}

# Flip strand flags and create an BAM index file
function flip_strands {
    samtools view -h ${BAM} ${REG} |
    awk 'BEGIN{ FS="\t"; OFS="\t"; }
    {
        if(substr($1,1,1) == "@") { # print BAM header
            print
        } else {
        if(and($2, 0x10)) { # flip "read reverse strand" flag
            $2 = $2-0x10
        } else {
        $2 = $2+0x10
    }
    print $_;
        }
    }' |
    samtools view -hb -o ${OUT}
    samtools index ${OUT}
}

# Print helpful usage information
function usage {
    >&2 echo -e "\nUsage: flip_pbat_strands.sh [-h,--help] [-r,--region] in_bam out_bam\n"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\tin_bam      : input BAM to flip strands"
    >&2 echo -e "\tout_bam     : output BAM to write flipped reads to\n"
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "\t-r,--region : region to flip, give as chr:start-end [default: ]"
    >&2 echo -e "\t-h,--help   : print help message and exit"
}

# Initialize default values for optional inputs
REG=""

# Process command line arguments
OPTS=$(getopt \
    --options hr: \
    --long help,region: \
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
        -r|--region )
            REG="${2}"
            shift 2
            ;;
        -- )
            shift
            break
            ;;
        * )
            >&2 echo "Unknown option: ${1}"
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

BAM=${1}
OUT=${2}

if [[ ! -f "${BAM}" ]]; then
    >&2 echo "Cannot find ${BAM}"
    exit 1
fi

if [[ -f ${OUT} ]] || [[ -f ${OUT}.bai ]]; then
    >&2 echo "${OUT} or ${OUT}.bai already exist"
    >&2 echo "Please delete and try again"
    exit 1
fi

flip_strands
