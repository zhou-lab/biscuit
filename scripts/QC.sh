#!/usr/bin/env bash
## make sure the following is in PATH
## biscuit samtools, bedtools, awk

function biscuitQC {

  mkdir -p $QCdir

  ## simple test to make sure all is available
  if [[ `which biscuit 2>&1 >/dev/null` ]]; then echo "biscuit does not exist in PATH"; exit 1; fi
  if [[ `which samtools 2>&1 >/dev/null` ]]; then echo "samtools does not exist in PATH"; exit 1; fi
  if [[ `which bedtools 2>&1 >/dev/null` ]]; then echo "bedtools does not exist in PATH"; exit 1; fi
  if [[ `which awk 2>&1 >/dev/null` ]]; then echo "awk does not exist in PATH"; exit 1; fi
  for var in BISCUIT_CPGBED BISCUIT_CGIBED BISCUIT_RMSK BISCUIT_EXON BISCUIT_GENE BISCUIT_TOPGC_BED BISCUIT_BOTGC_BED input_bam input_vcf; do
    if [[ ${!var} != "<unset>" && ! -f ${!var} ]]; then
      >&2 echo "$var: ${!var} does not exist."
      exit 1;
    fi
  done

  echo "Running."
  set -xe
  ##########################
  ## base coverage
  ##########################
  if [[ "$BISCUIT_QC_BASECOV" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_BASECOV ----"
    bedtools genomecov -bga -split -ibam $input_bam -g ${BISCUIT_REFERENCE}.fai | LC_ALL=C sort -k1,1 -k2,2n -T $QCdir >$QCdir/${sname}_bga.bed
    samtools view -q 40 -b $input_bam | bedtools genomecov -ibam stdin -g ${BISCUIT_REFERENCE}.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T $QCdir  >$QCdir/${sname}_bga_q40.bed
    awk '{cnt[$4]+=$3-$2}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_bga.bed | sort -k1,1n -T $QCdir >$QCdir/${sname}_bga_table.txt
    awk '{cnt[$4]+=$3-$2}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_bga_q40.bed | sort -k1,1n -T $QCdir >$QCdir/${sname}_bga_q40_table.txt
  fi
  
  ##########################
  ## duplicate_coverage
  ##########################
  [[ ! -f "$QCdir/${sname}_bga.bed" ]] && BISCUIT_QC_DUPLICATE=false
  [[ ! -f "$QCdir/${sname}_bga_q40.bed" ]] && BISCUIT_QC_DUPLCIATE=false
  if [[ "$BISCUIT_QC_DUPLICATE" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_DUPLICATE ----"
    # duplicate
    samtools view -f 0x400 -b $input_bam | bedtools genomecov -ibam stdin -g $BISCUIT_REFERENCE.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T $QCdir >$QCdir/${sname}_bga_dup.bed

    # duplication rate
    echo -ne "#bases covered by all reads: " >$QCdir/${sname}_dup_report.txt
    awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga.bed >>$QCdir/${sname}_dup_report.txt
    echo -ne "#bases covered by duplicate reads: " >>$QCdir/${sname}_dup_report.txt
    awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga_dup.bed >>$QCdir/${sname}_dup_report.txt

    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      # high GC content
      echo -ne "#high-GC bases covered by all reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga.bed -b $BISCUIT_TOPGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#high-GC bases covered by duplicate reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup.bed -b $BISCUIT_TOPGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt

      # low GC content
      echo -ne "#low-GC bases covered by all reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga.bed -b $BISCUIT_BOTGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#low-GC bases covered by duplicate reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup.bed -b $BISCUIT_BOTGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
    fi
    
    ## Q40
    # duplicate
    samtools view -f 0x400 -q 40 -b $input_bam | bedtools genomecov -ibam stdin -g $BISCUIT_REFERENCE.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T $QCdir >$QCdir/${sname}_bga_dup_q40.bed

    # duplication rate
    echo -ne "#bases covered by all q40-reads: " >>$QCdir/${sname}_dup_report.txt
    awk '$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga_q40.bed >>$QCdir/${sname}_dup_report.txt
    echo -ne "#bases covered by duplicate q40-reads: " >>$QCdir/${sname}_dup_report.txt
    awk '$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga_dup_q40.bed >>$QCdir/${sname}_dup_report.txt

    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      # high GC content
      echo -ne "#high-GC bases covered by all q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#high-GC bases covered by duplicate q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt

      # low GC content
      echo -ne "#low-GC bases covered by all q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#low-GC bases covered by duplicate q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
    fi
  fi
  
  ##########################
  ## cpg coverage
  ##########################

  [[ ! -f "$BISCUIT_CPGBED" ]] && BISCUIT_QC_CPGCOV=false
  [[ ! -f "$QCdir/${sname}_bga.bed" ]] && BISCUIT_QC_CPGCOV=false
  [[ ! -f "$QCdir/${sname}_bga_q40.bed" ]] && BISCUIT_QC_CPGCOV=false
  if [[ "$BISCUIT_QC_CPGCOV" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_CPGCOV ----"
    bedtools intersect -a $BISCUIT_CPGBED -b $QCdir/${sname}_bga.bed -wo -sorted | bedtools groupby -g 1-3 -c 7 -o min >$QCdir/${sname}_cpg.bed
    bedtools intersect -a $BISCUIT_CPGBED -b $QCdir/${sname}_bga_q40.bed -wo -sorted | bedtools groupby -g 1-3 -c 7 -o min >$QCdir/${sname}_cpg_q40.bed
    awk '{cnt[$4]+=1}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_cpg.bed | sort -k1,1n >$QCdir/${sname}_cpg_table.txt
    awk '{cnt[$4]+=1}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_cpg_q40.bed | sort -k1,1n >$QCdir/${sname}_cpg_q40_table.txt
  fi

  ##########################
  ## cpg distribution
  ##########################

  [[ ! -f "$QCdir/${sname}_cpg_q40.bed" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$QCdir/${sname}_cpg.bed" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_EXON" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_RMSK" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_GENE" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_CGIBED" ]] && BISCUIT_QC_CPGDIST=false
  if [[ "$BISCUIT_QC_CPGDIST" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_CPGDIST ----"
    # whole genome
    wc -l $QCdir/${sname}_cpg_q40.bed | awk -F" " '{printf("Territory\tAll\tUniqCov\tAllCov\nTotalCpGs\t%s",$1)}' >$QCdir/${sname}_cpg_dist_table.txt
    awk '$4>0{a+=1}END{printf("\t%d",a)}' $QCdir/${sname}_cpg_q40.bed >>$QCdir/${sname}_cpg_dist_table.txt
    awk '$4>0{a+=1}END{printf("\t%d\n",a)}' $QCdir/${sname}_cpg.bed >>$QCdir/${sname}_cpg_dist_table.txt

    # exon
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_EXON -sorted | wc -l | awk -F" " '{printf("ExonicCpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_EXON -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b $BISCUIT_EXON -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt

    # repeat
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_RMSK -sorted | wc -l | awk -F" " '{printf("RepeatCpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_RMSK -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b $BISCUIT_RMSK -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt

    # gene
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_GENE -sorted | wc -l | awk -F" " '{printf("GenicCpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_GENE -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b $BISCUIT_GENE -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt

    # CGI
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_CGIBED -sorted | wc -l | awk -F" " '{printf("CGICpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_CGIBED -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b $BISCUIT_CGIBED -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
  fi

  ##########################
  ## CGI coverage
  ##########################
  [[ ! -f "$BISCUIT_CGIBED" ]] && BISCUIT_QC_CGICOV=false
  [[ ! -f "$QCdir/${sname}_cpg_q40.bed" ]] && BISCUIT_QC_CGICOV=false
  if [[ "$BISCUIT_QC_CGICOV" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_CGICOV ----"
    # how CGI is covered by at least one q40-read in at least one CpG
    echo >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands\t" >>$QCdir/${sname}_cpg_dist_table.txt
    cat $BISCUIT_CGIBED | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least one CpG\t" >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_CGIBED -sorted -wo | awk '$4>0{print $5":"$6"-"$7}' | uniq -c | awk -F" " '{print $2"\t"$1}' | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least three CpGs\t" >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_CGIBED -sorted -wo | awk '$4>0{print $5":"$6"-"$7}' | uniq -c | awk -F" " '$1>=3{print $2"\t"$1}' | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least five CpGs\t" >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_CGIBED -sorted -wo | awk '$4>0{print $5":"$6"-"$7}' | uniq -c | awk -F" " '$1>=5{print $2"\t"$1}' | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least ten CpGs\t" >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_CGIBED -sorted -wo | awk '$4>0{print $5":"$6"-"$7}' | uniq -c | awk -F" " '$1>=10{print $2"\t"$1}' | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
  fi
  
  ##########################
  ## uniformity
  ##########################
  [[ ! -f "$QCdir/${sname}_bga.bed" ]] && BISCUIT_QC_UNIFORMITY=false
  [[ ! -f "$QCdir/${sname}_bga_q40.bed" ]] && BISCUIT_QC_UNIFORMITY=false
  if [[ "$BISCUIT_QC_UNIFORMITY" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_UNIFORMITY ----"
    awk -v sname="${sname}" '{cnt[$1]=$2}END{for (cov in cnt) {sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print "sample\tmu\tsigma\tcv\n"sname"_all\t"mu"\t"sigma"\t"sigma/mu}' $QCdir/${sname}_bga_q40_table.txt >$QCdir/${sname}_all_cv_table.txt
    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk -v sname="${sname}" -v output=$QCdir/${sname}_all_cv_table.txt '{cnt[$4]+=$3-$2}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_all_topgc\t"mu"\t"sigma"\t"sigma/mu >>output}' | sort -k1,1n >$QCdir/${sname}_coverage_cnts_topgc.txt
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk -v sname="${sname}" -v output=$QCdir/${sname}_all_cv_table.txt '{cnt[$4]+=$3-$2}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_all_botgc\t"mu"\t"sigma"\t"sigma/mu >>output}' | sort -k1,1n >$QCdir/${sname}_coverage_cnts_botgc.txt
    fi
  fi
  
  ##########################
  ## cpg uniformity
  ##########################
  [[ ! -f "$QCdir/${sname}_cpg_q40_table.txt" ]] && BISCUIT_QC_CPGUNIF=false
  [[ ! -f "$QCdir/${sname}_cpg_q40.bed" ]] && BISCUIT_QC_CPGUNIF=false
  if [[ "$BISCUIT_QC_CPGUNIF" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_CPGUNIF ----"
    awk -v sname="${sname}" '{cnt[$1]=$2}END{for(cov in cnt) {sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print "sample\tmu\tsigma\tcv\n"sname"_cpg\t"mu"\t"sigma"\t"sigma/mu}' $QCdir/${sname}_cpg_q40_table.txt >$QCdir/${sname}_cpg_cv_table.txt
    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk -v sname="${sname}" '{cnt[$4]+=1}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_cpg_topgc\t"mu"\t"sigma"\t"sigma/mu >>"$QCdir/${sname}_cpg_cv_table.txt"}' | sort -k1,1n >$QCdir/${sname}_cpg_coverage_cnts_topgc.txt
      bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk -v sname="${sname}" '{cnt[$4]+=1}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_cpg_botgc\t"mu"\t"sigma"\t"sigma/mu >>"$QCdir/${sname}_cpg_cv_table.txt"}' | sort -k1,1n >$QCdir/${sname}_cpg_coverage_cnts_botgc.txt
    fi
  fi
  
  ##########################
  ## bisulfite conversion
  ##########################
  [[ ! -f "$input_vcf" ]] && BISCUIT_QC_BSCONV=false
  if [[ "$BISCUIT_QC_BSCONV" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_BSCONV ----"
    samtools view -h -q 40 $input_bam | biscuit bsconv $BISCUIT_REFERENCE - | awk 'match($0,/ZN:Z:([^ ]*)/,a){print gensub(/[A-Z,_]+/, "\t", "g", a[1])}' | cut -f2,4,6,8 | awk -v OFS="\t" '{ra[$1]+=1;rc[$2]+=1;rg[$3]+=1;rt[$4]+=1;}END{for(k in ra) {print "CA", k, ra[k]} for(k in rc) {print "CC", k, rc[k]} for(k in rg) {print "CG", k, rg[k]} for(k in rt) {print "CT", k, rt[k]}}' | sort -k1,1 -k2,2n | awk 'BEGIN{print "CTXT\tnumRET\tCnt"}{print}' > $QCdir/${sname}_freqOfTotalRetentionPerRead.txt
    biscuit vcf2bed -et c $input_vcf | awk '{beta_sum[$6]+=$8; beta_cnt[$6]+=1;} END{print "CA\tCC\tCG\tCT"; print beta_sum["CA"]/beta_cnt["CA"]"\t"beta_sum["CC"]/beta_cnt["CC"]"\t"beta_sum["CG"]/beta_cnt["CG"]"\t"beta_sum["CT"]/beta_cnt["CT"];}' >$QCdir/${sname}_totalBaseConversionRate.txt
    samtools view -hq 40 $input_bam | biscuit bsconv -b $BISCUIT_REFERENCE - | awk '{for(i=1;i<=8;++i) a[i]+=$i;}END{print "CpA\tCpC\tCpG\tCpT"; print a[1]/(a[1]+a[2])"\t"a[3]/(a[3]+a[4])"\t"a[5]/(a[5]+a[6])"\t"a[7]/(a[7]+a[8]);}' >$QCdir/${sname}_totalReadConversionRate.txt
    samtools view -hq 40 $input_bam | biscuit cinread $BISCUIT_REFERENCE - -t ch -p QPAIR,CQPOS,CRETENTION | sort | uniq -c | awk -F" " '$4!="N"{print $2"\t"$3"\t"$4"\t"$1}' | sort -k1,1 -k2,2n -T $QCdir >$QCdir/${sname}_CpHRetentionByReadPos.txt
  fi

  ####################
  ## mapping_summary
  ####################
  if [[ "$BISCUIT_QC_MAPPING" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_MAPPING ----"
    biscuit cinread -p STRAND,BSSTRAND $BISCUIT_REFERENCE $input_bam | awk '{a[$1$2]+=1}END{for(strand in a) {print "strand\t"strand"\t"a[strand];}}' >$QCdir/${sname}_mapping_summary.txt
    samtools view -F 0x100 -f 0x4 $input_bam | wc -l | cat <(echo -ne "unmapped\t") - >$QCdir/${sname}_mapq_table.txt
    samtools view -F 0x104 $input_bam | awk '{cnt[$5]+=1}END{for(mapq in cnt) {print mapq"\t"cnt[mapq];}}' | sort -k1,1n >>$QCdir/${sname}_mapq_table.txt

  fi
}

input_vcf="<unset>"
QCdir="BISCUITqc"

function usage {
  >&2 echo "Usage: QC.sh [options] setup_file sample_name input_bam"
  >&2 echo "Input options:"
  >&2 echo "   -v|--vcf     vcf outupt from BISCUIT"
  >&2 echo "   -o|--output  output directory [$QCdir]"
  >&2 echo "   -h|--help    show this help"
  exit 1;
}

[[ "$#" -lt 3 ]] && usage;
input_bam="${@: -1}"
set -- "${@:1:$(($#-1))}"
sname="${@: -1}"
set -- "${@:1:$(($#-1))}"
setup_file="${@: -1}"
set -- "${@:1:$(($#-1))}"

if [[ ! -f "$setup_file" ]]; then
  echo "Setup file missing: $setup_file.";
  exit 1;
fi

source "$setup_file"

if [[ ! -f "${BISCUIT_REFERENCE}.fai" ]]; then
  >&2 echo "Cannot locate fai-indexed reference: ${BISCUIT_REFERENCE}.fai"
  >&2 echo "Please set resource links in the script.";
  exit 1;
fi

TEMP=`getopt -o v:o:h --long vcf:,sample::,help -n 'QC.sh' -- "$@"`
eval set -- "$TEMP"

while true ; do
  case "$1" in
    -v|--vcf) input_vcf=$2; shift 2;;
    -o|--output)
      case "$2" in
        "") shift 2;;
        *) QCdir=$2; shift 2;;
      esac;;
    -h|--help) usage; shift; break;;
    --) shift; break;;
    *) sample_name=$1; shift; input_bam=$2; break;;
  esac
done

>&2 echo "## Running BISCUIT QC script with following configuration ##"
>&2 echo "=============="
>&2 echo "sample name: $sname"
>&2 echo "input bam:   $input_bam"
>&2 echo "input vcf:   $input_vcf"
>&2 echo "output dir:  $QCdir"
>&2 echo "REFERENCE:   $BISCUIT_REFERENCE"
>&2 echo "CPGBED:      $BISCUIT_CPGBED"
>&2 echo "CGIBED:      $BISCUIT_CGIBED"
>&2 echo "RMSK:        $BISCUIT_RMSK"
>&2 echo "EXON:        $BISCUIT_EXON"
>&2 echo "GENE:        $BISCUIT_GENE"
>&2 echo "TOPGC_BED:   $BISCUIT_TOPGC_BED"
>&2 echo "BOTGC_BED:   $BISCUIT_BOTGC_BED"
>&2 echo "=============="
biscuitQC
>&2 echo -e "\nDone."
