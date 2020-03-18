#!/usr/bin/env bash
################################################################################
##
## Quality Control script for BISCUIT output
##
## Output from this script can be fed into MultiQC to produce a nice HTML output
## showing the different BISCUIT QC metrics
##
## Notes:
##   1.) biscuit, samtools, bedtools, and awk all must be in PATH for script to
##       work
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
##     - Bug fix in Coefficient of Variantion calculation
##
################################################################################

set -euo pipefail

# Check for biscuit, samtools, bedtools, awk in PATH
function check_path {
  if [[ `which biscuit 2>&1 > /dev/null` ]]; then
      >&2 echo "biscuit does not exist in PATH"
      exit 1
  fi
  if [[ `which samtools 2>&1 > /dev/null` ]]; then
      >&2 echo "samtools does not exist in PATH"
      exit 1
  fi
  if [[ `which bedtools 2>&1 > /dev/null` ]]; then
      >&2 echo "bedtools does not exist in PATH"
      exit 1
  fi
  if [[ `which awk 2>&1 > /dev/null` ]]; then
      >&2 echo "awk does not exist in PATH"
      exit 1
  fi
}

# Check that certain variables have been set and files exist
function check_variables {
    VARS="
    BISCUIT_CPGBED
    BISCUIT_CGIBED
    BISCUIT_RMSK
    BISCUIT_EXON
    BISCUIT_GENE
    BISCUIT_TOPGC_BED
    BISCUIT_BOTGC_BED
    in_vcf
    "

    for var in ${VARS}; do
        if [[ ${!var} != "<unused>" && ! -f ${!var} ]]; then
            >&2 echo "${var}: ${!var} does not exist"
            exit 1
        fi
    done
}

# Check if QC files have at least some information
function basic_check_output_filled {
    fileloc=${outdir}/${sample}
    TWO_LINE_FILES="
    ${fileloc}_all_cv_table.txt
    ${fileloc}_covdist_cpg_q40_botgc_table.txt
    ${fileloc}_covdist_cpg_q40_table.txt
    ${fileloc}_covdist_cpg_q40_topgc_table.txt
    ${fileloc}_covdist_cpg_table.txt
    ${fileloc}_covdist_q40_botgc_table.txt
    ${fileloc}_covdist_q40_table.txt
    ${fileloc}_covdist_q40_topgc_table.txt
    ${fileloc}_covdist_table.txt
    ${fileloc}_cpg_cv_table.txt
    ${fileloc}_cpg_dist_table.txt
    ${fileloc}_CpGRetentionByReadPos.txt
    ${fileloc}_CpGRetentionDist.txt
    ${fileloc}_CpHRetentionByReadPos.txt
    ${fileloc}_freqOfTotalRetentionPerRead.txt
    ${fileloc}_isize_score_table.txt
    ${fileloc}_mapq_table.txt
    ${fileloc}_totalBaseConversionRate.txt
    ${fileloc}_totalReadConversionRate.txt
    "
    ONE_LINE_FILES="
    ${fileloc}_dup_report.txt
    ${fileloc}_strand_table.txt
    "

    >&2 echo "Running basic check on BISCUIT QC output"
    >&2 echo "Will remove any files that were obviously not filled properly"
    >&2 echo "This is to avoid clashes when running MultiQC"

    # All files that have a description line, then a table header line
    for FILE in ${TWO_LINE_FILES}; do
        if [[ ! -f "${FILE}" ]]; then
            >&2 echo "----- ${FILE} ----- does not exist. Skipping!"
            continue
        fi
        if [[ `wc -l ${FILE} | awk '{print $1}'` -lt 3 ]]; then
            >&2 echo "----- ${FILE} ----- has no entries. Check file!"
            >&2 echo "Deleting ----- ${FILE} ----- to help with MultiQC report"
            rm -f ${FILE}
        fi
    done

    # Files with only a description line
    for FILE in ${ONE_LINE_FILES}; do
        if [[ ! -f "${FILE}" ]]; then
            >&2 echo "----- ${FILE} ----- was not initially created. Skipping!"
            continue
        fi
        if [[ `wc -l ${FILE} | awk '{print $1}'` -lt 2 ]]; then
            >&2 echo "----- ${FILE} ----- has no entries. Check file!"
            >&2 echo "Deleting ----- ${FILE} ----- to help with MultiQC report"
            rm -f ${FILE}
        fi
    done
}

# Workhorse function for processing BISCUIT QC
function biscuitQC {
    # TODO: Make this check an exit if it doesn't pass
    # TODO: Change output message to explain why it's not running
    # TODO?: Change to avoid occassional samtools view broken pipe error
    if [[ `samtools view ${in_bam} | head -n 5000 | wc -l` -le 4999 ]]; then
        >&2 echo "Less than 5000 reads, things are going to get sketchy"
    fi

    # Simple check for necessary command line tools
    check_path

    # Check variables and their associated files exist
    check_variables

    # Create outdir if it does not exist
    if [ ! -d ${outdir} ]; then
        mkdir -p ${outdir}
    fi

    >&2 echo -e "Starting BISCUIT QC at `date`\n"

    ###################################
    ## Base Coverage
    ###################################
    if [[ "${BISCUIT_QC_BASECOV}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_BASECOV ----- generated at `date`"

        # Find genome coverage of input BAM
        bedtools genomecov -bga -split -ibam ${in_bam} | \
        LC_ALL=C sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} \
            > ${outdir}/${sample}_genomecov_all.tmp.bed

        # Find genome coverage of input BAM with reads having MAPQ>40
        samtools view -q 40 -b ${in_bam} | \
        bedtools genomecov -bga -split -ibam stdin | \
        LC_ALL=C sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} \
            > ${outdir}/${sample}_genomecov_q40.tmp.bed

        echo -e "BISCUITqc Depth Distribution (All)" \
            > ${outdir}/${sample}_covdist_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_table.txt
        awk '{ cnt[$4] += $3-$2 } END {
            for (cov in cnt) { print int(cov)"\t"int(cnt[cov]) }
        }' ${outdir}/${sample}_genomecov_all.tmp.bed | \
        sort --parallel=${thread} -k1,1n -T ${outdir} \
            >> ${outdir}/${sample}_covdist_table.txt

        echo -e "BISCUITqc Depth Distribution (Q40)" \
            > ${outdir}/${sample}_covdist_q40_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_q40_table.txt
        awk '{ cnt[$4] += $3-$2 } END {
            for (cov in cnt) { print int(cov)"\t"int(cnt[cov]) }
        }' ${outdir}/${sample}_genomecov_q40.tmp.bed | \
        sort --parallel=${thread} -k1,1n -T ${outdir} \
            >> ${outdir}/${sample}_covdist_q40_table.txt
    fi

    #TODO: If all_dups or q40_dups are skipped, then fill in default values or
    #TODO: something so that those values aren't missing for MultiQC
    ###################################
    ## Duplicate Coverage
    ###################################
    if [[ ! -f "${outdir}/${sample}_genomecov_all.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_genomecov_all.tmp.bed"
        >&2 echo "Skipping duplicate QC"
        BISCUIT_QC_DUPLICATE=false
    elif [[ ! -f "${outdir}/${sample}_genomecov_q40.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_genomecov_q40.tmp.bed"
        >&2 echo "Skipping duplicate QC"
        BISCUIT_QC_DUPLICATE=false
    fi

    if [[ "${BISCUIT_QC_DUPLICATE}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_DUPLICATE ----- generated at `date`"

        # All duplicates
        # TODO?: Change to avoid occassional samtools view broken pipe error
        if [[ `samtools view -f 0x400 ${in_bam} | head -n 5 | wc -l` -ge 2 ]]; then
            samtools view -f 0x400 -b ${in_bam} |
            bedtools genomecov -bga -split -ibam stdin |
            LC_ALL=C sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} \
                > ${outdir}/${sample}_genomecov_all_dup.tmp.bed

            # Duplicate rate
            echo -e "BISCUITqc Read Duplication Table" \
                > ${outdir}/${sample}_dup_report.txt

            echo -ne "#bases covered by all reads: " \
                >> ${outdir}/${sample}_dup_report.txt
            awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                ${outdir}/${sample}_genomecov_all.tmp.bed \
                >> ${outdir}/${sample}_dup_report.txt
            echo -ne "#bases covered by duplicate reads: " \
                >> ${outdir}/${sample}_dup_report.txt
            awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                ${outdir}/${sample}_genomecov_all_dup.tmp.bed \
                >> ${outdir}/${sample}_dup_report.txt

            if [[ -f "${BISCUIT_TOPGC_BED}" && -f "${BISCUIT_BOTGC_BED}" ]]; then
                # High GC content
                echo -ne "#high-GC bases covered by all reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_all.tmp.bed \
                    -b ${BISCUIT_TOPGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt
                echo -ne "#high-GC bases covered by duplicate reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_all_dup.tmp.bed \
                    -b ${BISCUIT_TOPGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt

                # Low GC content
                echo -ne "#low-GC bases covered by all reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_all.tmp.bed \
                    -b ${BISCUIT_BOTGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt
                echo -ne "#low-GC bases covered by duplicate reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_all_dup.tmp.bed \
                    -b ${BISCUIT_BOTGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt
            fi
        fi

        # Q40 duplicates
        # TODO?: Change to avoid occassional samtools view broken pipe error
        if [[ `samtools view -f 0x400 -q 40 -b ${in_bam} | head -n 5 | wc -l` -ge 2 ]]; then
            samtools view -f 0x400 -q 40 -b ${in_bam} |
            bedtools genomecov -bga -split -ibam stdin |
            LC_ALL=C sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} \
                > ${outdir}/${sample}_genomecov_q40_dup.tmp.bed

            # Duplicate rate
            echo -ne "#bases covered by all q40-reads: " \
                >> ${outdir}/${sample}_dup_report.txt
            awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                ${outdir}/${sample}_genomecov_q40.tmp.bed \
                >> ${outdir}/${sample}_dup_report.txt
            echo -ne "#bases covered by duplicate q40-reads: " \
                >> ${outdir}/${sample}_dup_report.txt
            awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                ${outdir}/${sample}_genomecov_q40_dup.tmp.bed \
                >> ${outdir}/${sample}_dup_report.txt

            if [[ -f "${BISCUIT_TOPGC_BED}" && -f "${BISCUIT_BOTGC_BED}" ]]; then
                # High GC content
                echo -ne "#high-GC bases covered by all q40-reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_q40.tmp.bed \
                    -b ${BISCUIT_TOPGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt
                echo -ne "#high-GC bases covered by duplicate q40-reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_q40_dup.tmp.bed \
                    -b ${BISCUIT_TOPGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt

                # Low GC content
                echo -ne "#low-GC bases covered by all q40-reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_q40.tmp.bed \
                    -b ${BISCUIT_BOTGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt
                echo -ne "#low-GC bases covered by duplicate q40-reads: " \
                    >> ${outdir}/${sample}_dup_report.txt
                bedtools intersect -sorted \
                    -a ${outdir}/${sample}_genomecov_q40_dup.tmp.bed \
                    -b ${BISCUIT_BOTGC_BED} | \
                awk 'BEGIN { a=0 } $4>0 { a += $3-$2 } END { print a }' \
                    >> ${outdir}/${sample}_dup_report.txt
            fi
        fi
    fi

    ###################################
    ## CpG Coverage
    ###################################
    if [[ ! -f "${outdir}/${sample}_genomecov_all.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_genomecov_all.tmp.bed"
        >&2 echo "Skipping CpG coverage QC"
        BISCUIT_QC_CPGCOV=false
    elif [[ ! -f "${outdir}/${sample}_genomecov_q40.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_genomecov_q40.tmp.bed"
        >&2 echo "Skipping CpG coverage QC"
        BISCUIT_QC_CPGCOV=false
    fi

    if [[ -f "${BISCUIT_CPGBED}" && "${BISCUIT_QC_CPGCOV}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_CPGCOV ----- generated at `date`"

        # All CpGs with reads covering them
        bedtools intersect -sorted -wo \
            -a ${BISCUIT_CPGBED} \
            -b ${outdir}/${sample}_genomecov_all.tmp.bed | \
        bedtools groupby -g 1-3 -c 7 -o min \
            > ${outdir}/${sample}_cpg_all.tmp.bed

        # CpGs with reads having Q40 covering them
        bedtools intersect -sorted -wo \
            -a ${BISCUIT_CPGBED} \
            -b ${outdir}/${sample}_genomecov_q40.tmp.bed | \
        bedtools groupby -g 1-3 -c 7 -o min \
            > ${outdir}/${sample}_cpg_q40.tmp.bed

        echo -e "BISCUITqc CpG Depth Distribution (All)" \
            > ${outdir}/${sample}_covdist_cpg_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_cpg_table.txt
        awk '{ cnt[$4] += 1 } END {
            for (cov in cnt) { print int(cov)"\t"int(cnt[cov]) }
        }' ${outdir}/${sample}_cpg_all.tmp.bed | \
        sort --parallel=${thread} -k1,1n -T ${outdir} \
            >> ${outdir}/${sample}_covdist_cpg_table.txt

        echo -e "BISCUITqc CpG Depth Distribution (Q40)" \
            > ${outdir}/${sample}_covdist_cpg_q40_table.txt
        echo -e "depth\tcount" >> ${outdir}/${sample}_covdist_cpg_q40_table.txt
        awk '{ cnt[$4] += 1 } END {
            for (cov in cnt) { print int(cov)"\t"int(cnt[cov]) }
        }' ${outdir}/${sample}_cpg_q40.tmp.bed | \
        sort --parallel=${thread} -k1,1n -T ${outdir} \
            >> ${outdir}/${sample}_covdist_cpg_q40_table.txt
    fi

    ###################################
    ## CpG Distribution
    ###################################
    if [[ ! -f "${outdir}/${sample}_cpg_all.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_cpg_all.tmp.bed"
        >&2 echo "Skipping CpG coverage QC"
        BISCUIT_QC_CPGDIST=false
    elif [[ ! -f "${outdir}/${sample}_cpg_q40.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_cpg_q40.tmp.bed"
        >&2 echo "Skipping CpG coverage QC"
        BISCUIT_QC_CPGDIST=false
    fi

    if [[ "${BISCUIT_QC_CPGDIST}" == true ]] && \
       [[ -f "${BISCUIT_EXON}" ]] && \
       [[ -f "${BISCUIT_RMSK}" ]] && \
       [[ -f "${BISCUIT_GENE}" ]] && \
       [[ -f "${BISCUIT_CGIBED}" ]]; then
        >&2 echo "----- BISCUIT_QC_CPGDIST ----- generated at `date`"

        # Whole genome
        echo -e "BISCUITqc CpG Distribution Table" \
            > ${outdir}/${sample}_cpg_dist_table.txt
        echo -e "Territory\tAll\tUniqCov\tAllCov" \
            >> ${outdir}/${sample}_cpg_dist_table.txt

        wc -l ${outdir}/${sample}_cpg_q40.tmp.bed | \
        awk -F" " '{ printf("TotalCpGs\t%s", $1) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d", a) }' \
            ${outdir}/${sample}_cpg_q40.tmp.bed \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d\n", a) }' \
            ${outdir}/${sample}_cpg_all.tmp.bed \
            >> ${outdir}/${sample}_cpg_dist_table.txt

        # Exons
        bedtools merge -i ${BISCUIT_EXON} \
            > ${outdir}/${sample}_merged_exons.tmp.bed

        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_exons.tmp.bed | \
        wc -l | \
        awk -F" " '{ printf("ExonicCpGs\t%s", $1) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_exons.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_all.tmp.bed \
            -b ${outdir}/${sample}_merged_exons.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d\n", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt

        # Repeats
        bedtools merge -i ${BISCUIT_RMSK} \
            > ${outdir}/${sample}_merged_rmask.tmp.bed

        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_rmask.tmp.bed | \
        wc -l | \
        awk -F" " '{ printf("RepeatCpGs\t%s", $1) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_rmask.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_all.tmp.bed \
            -b ${outdir}/${sample}_merged_rmask.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d\n", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt

        # Gene
        bedtools merge -i ${BISCUIT_GENE} \
            > ${outdir}/${sample}_merged_genes.tmp.bed

        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_genes.tmp.bed | \
        wc -l | \
        awk -F" " '{ printf("GenicCpGs\t%s", $1) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_genes.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_all.tmp.bed \
            -b ${outdir}/${sample}_merged_genes.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d\n", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt

        # CpG Islands
        bedtools merge -i ${BISCUIT_CGIBED} \
            > ${outdir}/${sample}_merged_cgisl.tmp.bed

        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_cgisl.tmp.bed | \
        wc -l | \
        awk -F" " '{ printf("CGICpGs\t%s", $1) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_cgisl.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        bedtools intersect -sorted \
            -a ${outdir}/${sample}_cpg_all.tmp.bed \
            -b ${outdir}/${sample}_merged_cgisl.tmp.bed | \
        awk 'BEGIN { a=0 } $4>0 { a += 1 } END { printf("\t%d\n", a) }' \
            >> ${outdir}/${sample}_cpg_dist_table.txt

        # How many CGIs are covered by at least one q40-read in at least one CpG
        echo -ne "\n#CpG Islands\t" >> ${outdir}/${sample}_cpg_dist_table.txt
        zcat ${BISCUIT_CGIBED} | wc -l >> ${outdir}/${sample}_cpg_dist_table.txt

        bedtools intersect -sorted -wo \
            -a ${outdir}/${sample}_cpg_q40.tmp.bed \
            -b ${outdir}/${sample}_merged_cgisl.tmp.bed | \
        awk '$4>0 { print $5":"$6"-"$7 }' | \
        uniq -c | \
        awk -F" " '{ print $2"\t"$1 }' \
            >> ${outdir}/${sample}_cpg_dist_table.tmp.txt

        echo -ne "#CpG Islands covered by at least one q40-read in at least one CpG\t" \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        cat ${outdir}/${sample}_cpg_dist_table.tmp.txt | \
        wc -l >> ${outdir}/${sample}_cpg_dist_table.txt

        echo -ne "#CpG Islands covered by at least one q40-read in at least three CpGs\t" \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        awk -F" " '$2>=3'  ${outdir}/${sample}_cpg_dist_table.tmp.txt | \
        wc -l >> ${outdir}/${sample}_cpg_dist_table.txt

        echo -ne "#CpG Islands covered by at least one q40-read in at least five CpGs\t" \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        awk -F" " '$2>=5'  ${outdir}/${sample}_cpg_dist_table.tmp.txt | \
        wc -l >> ${outdir}/${sample}_cpg_dist_table.txt

        echo -ne "#CpG Islands covered by at least one q40-read in at least ten CpGs\t" \
            >> ${outdir}/${sample}_cpg_dist_table.txt
        awk -F" " '$2>=10'  ${outdir}/${sample}_cpg_dist_table.tmp.txt | \
        wc -l >> ${outdir}/${sample}_cpg_dist_table.txt
    fi

    ###################################
    ## Uniformity
    ###################################
    if [[ ! -f "${outdir}/${sample}_covdist_q40_table.txt" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_covdist_q40_table.txt"
        >&2 echo "Skipping unifomity QC"
        BISCUIT_QC_UNIFORMITY=false
    elif  [[ ! -f "${outdir}/${sample}_genomecov_q40.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_genomecov_q40.tmp.bed"
        >&2 echo "Skipping unifomity QC"
        BISCUIT_QC_UNIFORMITY=false
    fi

    if [[ "${BISCUIT_QC_UNIFORMITY}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_UNIFORMITY ----- generated at `date`"

        echo -e "BISCUITqc Uniformity Table" \
            > ${outdir}/${sample}_all_cv_table.txt
        echo -e "sample\tmu\tsigma\tcv" \
            >> ${outdir}/${sample}_all_cv_table.txt

        awk -v sname="${sample}" '{ cnt[$1] = $2 } END {
            for (cov in cnt) {
                sum_cov += cnt[cov]*cov
                sum_cnt += cnt[cov]
            }
            mu = sum_cov / sum_cnt
            for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) } 
            sigma = sqrt(sum_var / sum_cnt)
            print sname"_all\t"mu"\t"sigma"\t"sigma/mu
        }' ${outdir}/${sample}_covdist_q40_table.txt \
            >> ${outdir}/${sample}_all_cv_table.txt

        if [[ -f "${BISCUIT_TOPGC_BED}" && -f "${BISCUIT_BOTGC_BED}" ]]; then
            echo -e "BISCUITqc Depth Distribution (high GC, Q40)" \
                > ${outdir}/${sample}_covdist_q40_topgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_q40_topgc_table.txt

            bedtools intersect -sorted \
                -a ${outdir}/${sample}_genomecov_q40.tmp.bed \
                -b ${BISCUIT_TOPGC_BED} | \
            awk -v sname="${sample}" \
                -v output="${outdir}/${sample}_all_cv_table.txt" \
                '{ cnt[$4] += $3-$2 } END {
                    for (cov in cnt) {
                        print cov"\t"cnt[cov]
                        sum_cov += cnt[cov]*cov;
                        sum_cnt += cnt[cov]
                    }
                    mu = sum_cov / sum_cnt
                    for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) } 
                    sigma = sqrt(sum_var / sum_cnt)
                    print sname"_all_topgc\t"mu"\t"sigma"\t"sigma/mu >> output
                }' | \
            sort --parallel=${thread} -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_q40_topgc_table.txt

            echo -e "BISCUITqc Depth Distribution (low GC, Q40)" \
                > ${outdir}/${sample}_covdist_q40_botgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_q40_botgc_table.txt

            bedtools intersect -sorted \
                -a ${outdir}/${sample}_genomecov_q40.tmp.bed \
                -b ${BISCUIT_BOTGC_BED} | \
            awk -v sname="${sample}" \
                -v output="${outdir}/${sample}_all_cv_table.txt" \
                '{ cnt[$4] += $3-$2 } END {
                    for (cov in cnt) {
                        print cov"\t"cnt[cov]
                        sum_cov += cnt[cov]*cov;
                        sum_cnt += cnt[cov]
                    }
                    mu = sum_cov / sum_cnt
                    for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) } 
                    sigma = sqrt(sum_var / sum_cnt)
                    print sname"_all_botgc\t"mu"\t"sigma"\t"sigma/mu >> output
                }' | \
            sort --parallel=${thread} -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_q40_botgc_table.txt
        fi
    fi

    ###################################
    ## CpG Uniformity
    ###################################
    if [[ ! -f "${outdir}/${sample}_covdist_cpg_q40_table.txt" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_covdist_cpg_q40_table.txt"
        >&2 echo "Skipping CpG unifomity QC"
        BISCUIT_QC_CPGUNIF=false
    elif  [[ ! -f "${outdir}/${sample}_cpg_q40.tmp.bed" ]]; then
        >&2 echo "Missing ${outdir}/${sample}_cpg_q40.tmp.bed"
        >&2 echo "Skipping CpG unifomity QC"
        BISCUIT_QC_CPGUNIF=false
    fi

    if [[ "${BISCUIT_QC_CPGUNIF}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_CPGUNIF ----- generated at `date`"

        echo -e "BISCUITqc CpG Uniformity Table" \
            > ${outdir}/${sample}_cpg_cv_table.txt
        echo -e "sample\tmu\tsigma\tcv" \
            >> ${outdir}/${sample}_cpg_cv_table.txt

        awk -v sname="${sample}" '{ cnt[$1] = $2 } END {
            for (cov in cnt) {
                sum_cov += cnt[cov]*cov
                sum_cnt += cnt[cov]
            }
            mu = sum_cov / sum_cnt
            for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) } 
            sigma = sqrt(sum_var / sum_cnt)
            print sname"_cpg\t"mu"\t"sigma"\t"sigma/mu
        }' ${outdir}/${sample}_covdist_cpg_q40_table.txt \
            >> ${outdir}/${sample}_cpg_cv_table.txt

        if [[ -f "${BISCUIT_TOPGC_BED}" && -f "${BISCUIT_BOTGC_BED}" ]]; then
            echo -e "BISCUITqc CpG Depth Distribution (high GC, Q40)" \
                > ${outdir}/${sample}_covdist_cpg_q40_topgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_cpg_q40_topgc_table.txt

            bedtools intersect -sorted \
                -a ${outdir}/${sample}_cpg_q40.tmp.bed \
                -b ${BISCUIT_TOPGC_BED} | \
            awk -v sname="${sample}" \
                -v output="${outdir}/${sample}_cpg_cv_table.txt" \
                '{ cnt[$4] += 1 } END {
                    for (cov in cnt) {
                        print cov"\t"cnt[cov]
                        sum_cov += cnt[cov]*cov;
                        sum_cnt += cnt[cov]
                    }
                    mu = sum_cov / sum_cnt
                    for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) } 
                    sigma = sqrt(sum_var / sum_cnt)
                    print sname"_cpg_topgc\t"mu"\t"sigma"\t"sigma/mu >> output
                }' | \
            sort --parallel=${thread} -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_cpg_q40_topgc_table.txt

            echo -e "BISCUITqc CpG Depth Distribution (low GC, Q40)" \
                > ${outdir}/${sample}_covdist_cpg_q40_botgc_table.txt
            echo -e "depth\tcount" \
                >> ${outdir}/${sample}_covdist_cpg_q40_botgc_table.txt

            bedtools intersect -sorted \
                -a ${outdir}/${sample}_cpg_q40.tmp.bed \
                -b ${BISCUIT_BOTGC_BED} | \
            awk -v sname="${sample}" \
                -v output="${outdir}/${sample}_cpg_cv_table.txt" \
                '{ cnt[$4] += 1 } END {
                    for (cov in cnt) {
                        print cov"\t"cnt[cov]
                        sum_cov += cnt[cov]*cov;
                        sum_cnt += cnt[cov]
                    }
                    mu = sum_cov / sum_cnt
                    for (cov in cnt) { sum_var += cnt[cov]*((cov-mu)^2) } 
                    sigma = sqrt(sum_var / sum_cnt)
                    print sname"_cpg_botgc\t"mu"\t"sigma"\t"sigma/mu >> output
                }' | \
            sort --parallel=${thread} -k1,1n -T ${outdir} \
                >> ${outdir}/${sample}_covdist_cpg_q40_botgc_table.txt
        fi
    fi

    ###################################
    ## Bisulfite Conversion
    ###################################
    if [[ ! -f "${in_vcf}" ]]; then
        >&2 echo "Input VCF does not exist or was not supplied."
        >&2 echo "Skipping bisulfite conversion QC"
        BISCUIT_QC_BSCONV=false
    fi

    if [[ "${BISCUIT_QC_BSCONV}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_BSCONV ----- generated at `date`"

        echo -e "BISCUITqc Frequency of Total Retention per Read Table" \
            > ${outdir}/${sample}_freqOfTotalRetentionPerRead.txt

        samtools view -h -q 40 ${in_bam} | \
        biscuit bsconv ${genome} - ${outdir}/${sample}.tmp.bam &> /dev/null

        samtools view ${outdir}/${sample}.tmp.bam | \
        awk 'match($0, /ZN:Z:([^ ]*)/, a) {
            print gensub(/[A-Z,_]+/, "\t", "g", a[1])
        }' | \
        cut -f2,4,6,8 | \
        awk -v OFS="\t" \
            '{ ra[$1] += 1; rc[$2] += 1; rg[$3] += 1; rt[$4] += 1; } END {
                for (k in ra) { print "CA", k, ra[k] }
                for (k in rc) { print "CC", k, rc[k] }
                for (k in rg) { print "CG", k, rg[k] }
                for (k in rt) { print "CT", k, rt[k] }
        }' | \
        sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} | \
        awk 'BEGIN{ print "CTXT\tnumRET\tCnt" }{ print }' \
            >> ${outdir}/${sample}_freqOfTotalRetentionPerRead.txt

        echo -e "BISCUITqc Conversion Rate by Base Average Table" \
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
        }' \
            >> ${outdir}/${sample}_totalBaseConversionRate.txt

        echo -e "BISCUITqc Conversion Rate by Read Average Table" \
            > ${outdir}/${sample}_totalReadConversionRate.txt
        samtools view -h -q 40 -F 0x900 ${in_bam} |
        biscuit bsconv -b ${genome} - > ${outdir}/${sample}.tmp.txt

        cat ${outdir}/${sample}.tmp.txt | \
        awk '{ for (i=1;i<=8;++i) { a[i] += $i } } END {
            ca = a[1] / (a[1] + a[2])
            cc = a[3] / (a[3] + a[4])
            cg = a[5] / (a[5] + a[6])
            ct = a[7] / (a[7] + a[8])
            print "CpA\tCpC\tCpG\tCpT"
            print ca"\t"cc"\t"cg"\t"ct
        }' >> ${outdir}/${sample}_totalReadConversionRate.txt

        echo -e "BISCUITqc CpH Retention by Read Position Table" \
            > ${outdir}/${sample}_CpHRetentionByReadPos.txt
        echo -e "ReadInPair\tPosition\tConversion/Retention\tCount" \
            >> ${outdir}/${sample}_CpHRetentionByReadPos.txt
        samtools view -h -q 40 ${in_bam} | \
        biscuit cinread ${genome} - -t ch -p QPAIR,CQPOS,CRETENTION | \
        sort --parallel=${thread} -T ${outdir} | \
        uniq -c | \
        awk -F" " '$4 != "N" { print $2"\t"$3"\t"$4"\t"$1 }' | \
        sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} \
            >> ${outdir}/${sample}_CpHRetentionByReadPos.txt

        echo -e "BISCUITqc CpG Retention by Read Position Table" \
            > ${outdir}/${sample}_CpGRetentionByReadPos.txt
        echo -e "ReadInPair\tPosition\tConversion/Retention\tCount" \
            >> ${outdir}/${sample}_CpGRetentionByReadPos.txt
        samtools view -h -q 40 ${in_bam} | \
        biscuit cinread ${genome} - -t cg -p QPAIR,CQPOS,CRETENTION | \
        sort --parallel=${thread} -T ${outdir} | \
        uniq -c | \
        awk -F" " '$4 != "N" { print $2"\t"$3"\t"$4"\t"$1 }' | \
        sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} \
            >> ${outdir}/${sample}_CpGRetentionByReadPos.txt
    fi

    ###################################
    ## Mapping Summary
    ###################################
    if [[ "${BISCUIT_QC_MAPPING}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_MAPPING ----- generated at `date`"

        echo -e "BISCUITqc Strand Table" > ${outdir}/${sample}_strand_table.txt
        biscuit cinread ${genome} ${in_bam} -p QPAIR,STRAND,BSSTRAND | \
        awk '{ a[$1$2$3] += 1 } END {
            for (strand in a) { print "strand\t"strand"\t"a[strand] } }' \
            >> ${outdir}/${sample}_strand_table.txt

        echo -e "BISCUITqc Mapping Quality Table" \
            > ${outdir}/${sample}_mapq_table.txt
        echo -e "MapQ\tCount" >> ${outdir}/${sample}_mapq_table.txt
        samtools view -F 0x100 -f 0x4 ${in_bam} | \
        wc -l | \
        cat <(echo -ne "unmapped\t") - >> ${outdir}/${sample}_mapq_table.txt

        samtools view -F 0x104 ${in_bam} | \
        awk '{ cnt[$5] += 1 } END {
            for (mapq in cnt) { print mapq"\t"cnt[mapq] } }' | \
        sort --parallel=${thread} -k1,1n -T ${outdir} \
            >> ${outdir}/${sample}_mapq_table.txt

        # Insert size - excludes reads by AS (40) and MAPQ (40)
        echo -e "BISCUITqc Insert Size, Score Table" \
            > ${outdir}/${sample}_isize_score_table.txt
        echo -e "InsertSize/Score\tValue\tFraction" \
            >> ${outdir}/${sample}_isize_score_table.txt
        samtools view -F 0x104 ${in_bam} | \
        awk '{ match($0, /AS:i:([0-9]*)/, a); score[a[1]] += 1; sumscore += 1
            if (and($2, 0x2) && a[1] >= 40 && $5 >= 40 && $9 >= 0 && $9 <= 2000) {
                isize[$9] += 1; sumisize += 1
            } } END {
                for (k in isize) { print "I\t"k"\t"isize[k]/sumisize }
                for (k in score) { print "S\t"k"\t"score[k]/sumscore }
        }' | \
        sort --parallel=${thread} -k1,1 -k2,2n -T ${outdir} \
            >> ${outdir}/${sample}_isize_score_table.txt
    fi

    ###################################
    ## CpG Retention Distribution
    ###################################
    if [[ ! -f "${in_vcf}" ]]; then
        >&2 echo "Input VCF does not exist or was not supplied."
        >&2 echo "Skipping CpG retention distribution QC"
        BISCUIT_QC_BETAS=false
    fi

    if [[ "${BISCUIT_QC_BETAS}" == true ]]; then
        >&2 echo "----- BISCUIT_QC_BETAS ----- generated at `date`"

        echo -e "BISCUITqc Retention Distribution Table" \
            > ${outdir}/${sample}_CpGRetentionDist.txt
        echo -e "RetentionFraction\tCount" \
            >> ${outdir}/${sample}_CpGRetentionDist.txt
        biscuit vcf2bed -k 1 -t cg ${in_vcf} | \
        awk '{
            if ($5>=3) { a[sprintf("%3.0f", $4*100)] += 1 }
            else { low_cov += 1 } } END {
                print "-1\t"low_cov
                for (beta in a) { print beta"\t"a[beta] }
        }' | \
        sort --parallel=${thread} -k1,1n -T ${outdir} \
            >> ${outdir}/${sample}_CpGRetentionDist.txt
    fi


    ###################################
    ## Remove temporary files
    ###################################
    if [[ ! "${keep_tmp}" == true ]]; then
        rm -f ${outdir}/${sample}*.tmp.*
    fi

    ###################################
    ## Run basic checks on output files
    ###################################
    basic_check_output_filled 

    >&2 echo -e "\nFinished BISCUIT QC at `date`"
}


# Print helpful usage information
function usage {
    >&2 echo -e "\nUsage: QC.sh [-h,--help] [-v,--vcf] [-o,--outdir] [-t,--threads] [-k,--keep-tmp-files] assets_directory genome sample_name in_bam\n"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\tassets_directory    : Path to assets directory"
    >&2 echo -e "\tgenome              : Path to reference FASTA file used in alignment"
    >&2 echo -e "\tsample_name         : Prefix of output files"
    >&2 echo -e "\tinput_bam           : Aligned BAM from BISCUIT\n"
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "\t-h,--help           : Print help message and exit"
    >&2 echo -e "\t-v,--vcf            : Path to VCF output from BISCUIT [DEFAULT: <unused>]"
    >&2 echo -e "\t-o,--outdir         : Output directory [DEFAULT: BISCUITqc]"
    >&2 echo -e "\t-t,--threads        : Number of threads to use [DEFAULT: 1]"
    >&2 echo -e "\t-k,--keep-tmp-files : Flag to keep temporary files for debugging [DEFAULT: Remove files]\n"
}

# Initialize default values for optional inputs
in_vcf="<unused>"
outdir="BISCUITqc"
thread=1
keep_tmp=false

# Process command line arguments
OPTS=$(getopt \
    --options hv:o:t:k \
    --long help,vcf:,outdir:,threads:,keep-bed-files \
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
        -v|--vcf )
            in_vcf="$2"
            shift 2
            ;;
        -o|--outdir )
            outdir="$2"
            shift 2
            ;;
        -t|--threads )
            thread="$2"
            shift 2
            ;;
        -k|--keep-tmp-files )
            keep_tmp=true
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

source $(dirname ${BASH_SOURCE[0]})/setup.sh $assets

>&2 echo "## Running BISCUIT QC script with following configuration ##"
>&2 echo "=============="
>&2 echo "Sample Name        : ${sample}"
>&2 echo "Input BAM          : ${in_bam}"
>&2 echo "Input VCF          : ${in_vcf}"
>&2 echo "Output Directory   : ${outdir}"
>&2 echo "# sort threads     : ${thread}"
>&2 echo "Assets Directory   : ${assets}"
>&2 echo "Reference          : ${genome}"
>&2 echo "Keep *.tmp.* files : ${keep_tmp}"
>&2 echo "CPGBED             : ${BISCUIT_CPGBED}"
>&2 echo "CGIBED             : ${BISCUIT_CGIBED}"
>&2 echo "RMSK               : ${BISCUIT_RMSK}"
>&2 echo "EXON               : ${BISCUIT_EXON}"
>&2 echo "GENE               : ${BISCUIT_GENE}"
>&2 echo "TOPGC_BED          : ${BISCUIT_TOPGC_BED}"
>&2 echo "BOTGC_BED          : ${BISCUIT_BOTGC_BED}"
>&2 echo "=============="
biscuitQC
