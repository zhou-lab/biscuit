#!/usr/bin/env bash
################################################################################
##
## Quality Control script for BISCUIT output
##
## Output from this script can be fed into MultiQC to produce a nice HTML output
## showing the different BISCUIT QC metrics
##
## Notes:
##   1.) biscuit, samtools, bedtools, and awk all must be in PATH for
##       script to work
##
## Created by:
##   Wanding Zhou
##
## Creation date:
##   May 2019
##
## Update notes:
##   Dec 2019 -
##     - Clean up code to make more readable
##     - Catch empty files, alert user, and remove files
##   Mar 2020 -
##     - Reworking of QC script
##         -- Removed repeats of functions and code
##         -- Removed python references for using BASH throughout
##         -- All sorts now have parallel option
##         -- Consistent tmp file naming scheme
##         -- <unset> is now <unused>
##         -- Assorted other changes
##     - Bug fix in Coefficient of Variation calculation
##   Apr 2020 -
##     - Refactoring QC script to reduce time spent running
##   May 2020 -
##     - Changing flags as necessary for updates to subcommand flags
##     - Adding PATH check for GNU parallel
##   Jun 2020 -
##     - Adding PATH check for GNU awk
##   Jul 2021 -
##     - Create qc subcommand and integrate into QC.sh
##   Jan 2022 -
##     - Remove GNU parallel dependency
##
################################################################################

set -euo pipefail

# Check for biscuit, samtools, bedtools, awk in PATH
function check_path {
  if [[ `which biscuit 2>&1 > /dev/null` ]]; then
      >&2 echo "biscuit does not exist in PATH"
      exit 1
  else
      >&2 echo "Using biscuit found at: `which biscuit`"
  fi
  if [[ `which samtools 2>&1 > /dev/null` ]]; then
      >&2 echo "samtools does not exist in PATH"
      exit 1
  else
      >&2 echo "Using samtools found at: `which samtools`"
  fi
  if [[ `which bedtools 2>&1 > /dev/null` ]]; then
      >&2 echo "bedtools does not exist in PATH"
      exit 1
  else
      >&2 echo "Using bedtools found at: `which bedtools`"
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

# Check that certain variables have been set and files exist
function check_variables {
    VARS="
    BISCUIT_CPGS
    BISCUIT_TOPGC
    BISCUIT_BOTGC
    in_vcf
    "

    for var in ${VARS}; do
        if [[ ${!var} != "<unused>" && ! -f ${!var} ]]; then
            >&2 echo "${var}: ${!var} does not exist"
            exit 1
        fi
    done
}

# Wait for jobs, exit if they're unsuccessful
# Via: https://stackoverflow.com/questions/1131484/wait-for-bash-background-jobs-in-script-to-be-finished/
function wait_for_jobs {
    while true; do
	wait -n || {
	    code="$?"
	    ([[ $code = "127" ]] && exit 0 || exit "$code")
	    break
	}
    done;
}

# Workhorse function for processing BISCUIT QC
function biscuitQC {
    # Simple check for necessary command line tools
    check_path

    # Check variables and their associated files exist
    #check_variables

    # Create outdir if it does not exist
    if [ ! -d ${outdir} ]; then
        mkdir -p ${outdir}
    fi

    >&2 echo -e "Starting BISCUIT QC at `date`\n"

    # MAPQ, Insert Size, Duplicate Rates, Strand Info, and Read-averaged conversion
    if [[ "${single_end}" == true ]]; then
        biscuit qc -s ${genome} ${in_bam} ${outdir}/${sample}
    else
        biscuit qc ${genome} ${in_bam} ${outdir}/${sample}
    fi

    if [[ "${run_cov_qc}" == true ]]; then
        # Create genomecov_all, genomecov_q40, genomecov_all_dup, genomecov_q40_dup
        # Spawn these to the background
        bedtools genomecov -bga -split -ibam ${in_bam} | LC_ALL=C sort -k1,1 -k2,2n -T ${outdir} > ${outdir}/${sample}_genomecov_all.tmp.bed &
        samtools view -q 40 -b ${in_bam} | bedtools genomecov -bga -split -ibam stdin | LC_ALL=C sort -k1,1 -k2,2n -T ${outdir} > ${outdir}/${sample}_genomecov_q40.tmp.bed &
        samtools view -f 0x400 -b ${in_bam} | bedtools genomecov -bga -split -ibam stdin | LC_ALL=C sort -k1,1 -k2,2n -T ${outdir} > ${outdir}/${sample}_genomecov_all_dup.tmp.bed &
        samtools view -f 0x400 -q 40 -b ${in_bam} | bedtools genomecov -bga -split -ibam stdin | LC_ALL=C sort -k1,1 -k2,2n -T ${outdir} > ${outdir}/${sample}_genomecov_q40_dup.tmp.bed &

        wait_for_jobs

        # Create cpg_all, cpg_q40
        bedtools intersect -sorted -wo -b ${outdir}/${sample}_genomecov_all.tmp.bed -a ${BISCUIT_CPGS} | bedtools groupby -g 1-3 -c 7 -o min > ${outdir}/${sample}_cpg_all.tmp.bed &
        bedtools intersect -sorted -wo -b ${outdir}/${sample}_genomecov_q40.tmp.bed -a ${BISCUIT_CPGS} | bedtools groupby -g 1-3 -c 7 -o min > ${outdir}/${sample}_cpg_q40.tmp.bed &

        wait_for_jobs

        # Coverage distributions and uniformity
        echo -e "BISCUITqc Uniformity Table" > ${outdir}/${sample}_cv_table.txt
        echo -e "group\tmu\tsigma\tcv" >> ${outdir}/${sample}_cv_table.txt

        echo -e "BISCUITqc Depth Distribution - All Bases" \
            > ${outdir}/${sample}_covdist_all_base_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_all_base_table.txt
        cat ${outdir}/${sample}_genomecov_all.tmp.bed | \
        awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += $3-$2 } END {
            for (cov in cnt) {
                print int(cov)"\t"int(cnt[cov])
                sum_cov += cnt[cov]*cov
                sum_cnt += cnt[cov]
            }
            mu = sum_cov / sum_cnt
            for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
            sigma = sqrt(sum_var / sum_cnt)
            print "all_base\t"mu"\t"sigma"\t"sigma/mu >> output
        }' | \
        sort -k1,1n -T ${outdir} >> ${outdir}/${sample}_covdist_all_base_table.txt

        echo -e "BISCUITqc Depth Distribution - All CpGs" \
            > ${outdir}/${sample}_covdist_all_cpg_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_all_cpg_table.txt
        cat ${outdir}/${sample}_cpg_all.tmp.bed | \
        awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += 1 } END {
            for (cov in cnt) {
                print int(cov)"\t"int(cnt[cov])
                sum_cov += cnt[cov]*cov
                sum_cnt += cnt[cov]
            }
            mu = sum_cov / sum_cnt
            for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
            sigma = sqrt(sum_var / sum_cnt)
            print "all_cpg\t"mu"\t"sigma"\t"sigma/mu >> output
        }' | \
        sort -k1,1n -T ${outdir} >> ${outdir}/${sample}_covdist_all_cpg_table.txt

        echo -e "BISCUITqc Depth Distribution - Q40 Bases" \
            > ${outdir}/${sample}_covdist_q40_base_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_q40_base_table.txt
        cat ${outdir}/${sample}_genomecov_q40.tmp.bed | \
        awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += $3-$2 } END {
            for (cov in cnt) {
                print int(cov)"\t"int(cnt[cov])
                sum_cov += cnt[cov]*cov
                sum_cnt += cnt[cov]
            }
            mu = sum_cov / sum_cnt
            for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
            sigma = sqrt(sum_var / sum_cnt)
            print "q40_base\t"mu"\t"sigma"\t"sigma/mu >> output
        }' | \
        sort -k1,1n -T ${outdir} >> ${outdir}/${sample}_covdist_q40_base_table.txt

        echo -e "BISCUITqc Depth Distribution - Q40 CpGs" \
            > ${outdir}/${sample}_covdist_q40_cpg_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_q40_cpg_table.txt
        cat ${outdir}/${sample}_cpg_q40.tmp.bed | \
        awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += 1 } END {
            for (cov in cnt) {
                print int(cov)"\t"int(cnt[cov])
                sum_cov += cnt[cov]*cov
                sum_cnt += cnt[cov]
            }
            mu = sum_cov / sum_cnt
            for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
            sigma = sqrt(sum_var / sum_cnt)
            print "q40_cpg\t"mu"\t"sigma"\t"sigma/mu >> output
        }' | \
        sort -k1,1n -T ${outdir} >> ${outdir}/${sample}_covdist_q40_cpg_table.txt

        if [[ -f "${BISCUIT_TOPGC}" && -f "${BISCUIT_BOTGC}" ]]; then
            echo -e "BISCUITqc Depth Distribution - All Top GC Bases" \
                > ${outdir}/${sample}_covdist_all_base_topgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_all_base_topgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_genomecov_all.tmp.bed \
                -b ${BISCUIT_TOPGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += $3-$2 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "all_base_topgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_all_base_topgc_table.txt

            echo -e "BISCUITqc Depth Distribution - All Top GC CpGs" \
                > ${outdir}/${sample}_covdist_all_cpg_topgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_all_cpg_topgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_cpg_all.tmp.bed \
                -b ${BISCUIT_TOPGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += 1 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "all_cpg_topgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_all_cpg_topgc_table.txt

            echo -e "BISCUITqc Depth Distribution - Q40 Top GC Bases" \
                > ${outdir}/${sample}_covdist_q40_base_topgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_q40_base_topgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_genomecov_q40.tmp.bed \
                -b ${BISCUIT_TOPGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += $3-$2 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "q40_base_topgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_q40_base_topgc_table.txt

            echo -e "BISCUITqc Depth Distribution - Q40 Top GC CpGs" \
                > ${outdir}/${sample}_covdist_q40_cpg_topgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_q40_cpg_topgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_cpg_q40.tmp.bed \
                -b ${BISCUIT_TOPGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += 1 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "q40_cpg_topgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_q40_cpg_topgc_table.txt

            echo -e "BISCUITqc Depth Distribution - All Bot GC Bases" \
                > ${outdir}/${sample}_covdist_all_base_botgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_all_base_botgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_genomecov_all.tmp.bed \
                -b ${BISCUIT_BOTGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += $3-$2 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "all_base_botgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_all_base_botgc_table.txt

            echo -e "BISCUITqc Depth Distribution - All Bot GC CpGs" \
                > ${outdir}/${sample}_covdist_all_cpg_botgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_all_cpg_botgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_cpg_all.tmp.bed \
                -b ${BISCUIT_BOTGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += 1 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "all_cpg_botgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_all_cpg_botgc_table.txt

            echo -e "BISCUITqc Depth Distribution - Q40 Bot GC Bases" \
                > ${outdir}/${sample}_covdist_q40_base_botgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_q40_base_botgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_genomecov_q40.tmp.bed \
                -b ${BISCUIT_BOTGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += $3-$2 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "q40_base_botgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_q40_base_botgc_table.txt

            echo -e "BISCUITqc Depth Distribution - Q40 Bot GC CpGs" \
                > ${outdir}/${sample}_covdist_q40_cpg_botgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_q40_cpg_botgc_table.txt
            bedtools intersect -sorted \
                -a ${outdir}/${sample}_cpg_q40.tmp.bed \
                -b ${BISCUIT_BOTGC} | \
            awk -v output="${outdir}/${sample}_cv_table.txt" '{ cnt[$4] += 1 } END {
                for (cov in cnt) {
                    print int(cov)"\t"int(cnt[cov])
                    sum_cov += cnt[cov]*cov
                    sum_cnt += cnt[cov]
                }
                mu = sum_cov / sum_cnt
                for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) }
                sigma = sqrt(sum_var / sum_cnt)
                print "q40_cpg_botgc\t"mu"\t"sigma"\t"sigma/mu >> output
            }' | \
            sort -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_q40_cpg_botgc_table.txt
        else
            >&2 echo -ne "Either ${BISCUIT_TOPGC} or ${BISCUIT_BOTGC} could not be "
            >&2 echo -ne "found. covdist files and uniformity metrics related to "
            >&2 echo -ne "top/bottom GC-content deciles will not be generated."
        fi
    fi

    if [[ -f ${in_vcf} ]]; then
        echo "BISCUITqc Conversion Rate by Base Average Table" \
            > ${outdir}/${sample}_totalBaseConversionRate.txt
        biscuit vcf2bed -e -t c ${in_vcf} | \
        awk '{ beta_sum[$6] += $8; beta_cnt[$6] += 1 } END {
            print "CA\tCC\tCG\tCT"
            if ( beta_cnt["CA"] < 20 ) {
                ca_frac = -1
            } else {
                ca_frac = beta_sum["CA"] / beta_cnt["CA"]
            }
            if ( beta_cnt["CC"] < 20 ) {
                cc_frac = -1
            } else {
                cc_frac = beta_sum["CC"] / beta_cnt["CC"]
            }
            if ( beta_cnt["CG"] < 20 ) {
                cg_frac = -1
            } else {
                cg_frac = beta_sum["CG"] / beta_cnt["CG"]
            }
            if ( beta_cnt["CT"] < 20 ) {
                ct_frac = -1
            } else {
                ct_frac = beta_sum["CT"] / beta_cnt["CT"]
            }
            print ca_frac"\t"cc_frac"\t"cg_frac"\t"ct_frac
        }' >> ${outdir}/${sample}_totalBaseConversionRate.txt
    else
        >&2 echo -ne "Input VCF not supplied. "
        >&2 echo -ne "${sample}_totalBaseConversionRate.txt will not be generated."
    fi

    ###################################
    ## Remove temporary files
    ###################################
    if [[ ! "${keep_tmp}" == true ]]; then
        rm -f ${outdir}/${sample}*.tmp.*
    fi

    >&2 echo -e "\nFinished BISCUIT QC at `date`"
}

# Print helpful usage information
function usage {
    >&2 echo -e "\nUsage: QC.sh [-h,--help] [-s,--single-end] [-v,--vcf] [-o,--outdir] [-k,--keep-tmp-files] [-n,--no-cov-qc] assets_directory genome sample_name in_bam\n"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\tassets_directory    : Path to assets directory"
    >&2 echo -e "\tgenome              : Path to reference FASTA file used in alignment"
    >&2 echo -e "\tsample_name         : Prefix of output files"
    >&2 echo -e "\tinput_bam           : Aligned BAM from BISCUIT\n"
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "\t-h,--help           : Print help message and exit"
    >&2 echo -e "\t-s,--single-end     : Input BAM is from single end data [DEFAULT: Assumes paired-end]"
    >&2 echo -e "\t-v,--vcf            : Path to VCF output from BISCUIT [DEFAULT: <unused>]"
    >&2 echo -e "\t-o,--outdir         : Output directory [DEFAULT: BISCUITqc]"
    >&2 echo -e "\t-k,--keep-tmp-files : Flag to keep temporary files for debugging [DEFAULT: Delete files]"
    >&2 echo -e "\t-n,--no-cov-qc      : Do not perform coverage or coverage uniformity QC [DEFAULT: Runs coverage QC]\n"
}

# Initialize default values for optional inputs
in_vcf="<unused>"
outdir="BISCUITqc"
keep_tmp=false
single_end=false
run_cov_qc=true

# Process command line arguments
OPTS=$(getopt \
    --options hsv:o:kn \
    --long help,single-end,vcf:,outdir:,keep-bed-files,no-cov-qc \
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
        -s|--single-end )
            single_end=true
            shift
            ;;
        -v|--vcf )
            in_vcf="$2"
            shift 2
            ;;
        -o|--outdir )
            outdir="$2"
            shift 2
            ;;
        -k|--keep-tmp-files )
            keep_tmp=true
            shift
            ;;
        -n|--no-cov-qc )
            run_cov_qc=false
            shift
            ;;
        -- )
            shift
            break
            ;;
        * )
            >&2 echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Make sure there are the correct number of inputs
if [[ $# -ne 4 ]]; then
    >&2 echo "$0: Missing inputs"
    usage
    exit 1
fi

# Fill required positional arguments
assets=$1
genome=$2
sample=$3
in_bam=$4

# Do some checks on the given inputs
if [[ ! -d "$assets" ]]; then
    >&2 echo "Assets directory missing: $assets"
    exit 1
fi

if [[ ! -f "${genome}.fai" ]]; then
    >&2 echo "Cannot locate fai-indexed reference: ${genome}.fai"
    >&2 echo "Please provide a viable path to the reference genome FASTA file."
    exit 1
fi

if [[ ! -f "${in_bam}" ]]; then
    >&2 echo "Cannot locate aligned BAM: ${in_bam}"
    >&2 echo "Please provide an existing aligned BAM."
    exit 1
fi

# Set variables for supplementary BED files
BISCUIT_CPGS="${assets}/cpg.bed.gz"
BISCUIT_TOPGC="${assets}/windows100bp.gc_content.top10p.bed.gz"
BISCUIT_BOTGC="${assets}/windows100bp.gc_content.bot10p.bed.gz"

>&2 echo "## Running BISCUIT QC script with following configuration ##"
>&2 echo "=============="
>&2 echo "Sample Name        : ${sample}"
>&2 echo "Input BAM          : ${in_bam}"
>&2 echo "Input VCF          : ${in_vcf}"
>&2 echo "Output Directory   : ${outdir}"
>&2 echo "Assets Directory   : ${assets}"
>&2 echo "Reference          : ${genome}"
>&2 echo "Keep *.tmp.* files : ${keep_tmp}"
>&2 echo "Run coverage QC    : ${run_cov_qc}"
>&2 echo "Single-end data    : ${single_end}"
>&2 echo "CPGS               : ${BISCUIT_CPGS}"
>&2 echo "TOPGC              : ${BISCUIT_TOPGC}"
>&2 echo "BOTGC              : ${BISCUIT_BOTGC}"
>&2 echo "=============="

biscuitQC
