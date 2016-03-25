
# bin/biscuit align ~/reference/mm10/mm10.fa /data/largeS2/pl-bs/2015-06-09-mouse-WGBS/run3_H7CVYBGXX/data/fastq/PL4-10-15WGBS1_L000_R1_001.fastq.gz /data/largeS2/pl-bs/2015-06-09-mouse-WGBS/run3_H7CVYBGXX/data/fastq/PL4-10-15WGBS1_L000_R2_001.fastq.gz -t 70 | samtools view -bS - > /data/largeS2/pl-bs/H7CVYBGXX_WGBS1.bam

export PROG=~/tools/biscuit/development/biscuit/bin/biscuit
export GENOME_DIR=~/references/

function pecho() {
  echo "[$(date)] Running:  "$@ >&2
}

function decho() {
  pecho $@
  eval $@
}

function biscuittest_pileup {
  base=test/Smadh3_chr19_chrM
  [[ -d $base/out_pileup ]] || mkdir $base/out_pileup
  rm -f $base/out_pileup/*
  decho "$PROG pileup -i $base/raw_bam/WGBS_Smadh3_chr19_chrM.bam -r $GENOME_DIR/mm10/mm10.fa -o $base/out_pileup/out.vcf -q 18"
  decho "bgzip $base/out_pileup/out.vcf"
  decho "tabix -p vcf $base/out_pileup/out.vcf.gz"

  base=test/NIC1254A46
  decho "$PROG pileup -r $GENOME_DIR/hg19/hg19.fa -i $base/raw_bam/NIC1254A46.bam -q 1 -g chr20:29570686-29570686 -v 3"
  decho "$PROG pileup -r $GENOME_DIR/hg19/hg19.fa -i $base/raw_bam/NIC1254A46.bam -q 1 -g chr20:26138808-26138808 -v 1 -n 999"
  decho "$PROG pileup -r $GENOME_DIR/hg19/hg19.fa -i $base/raw_bam/NIC1254A46.bam -q 1 -g chr20:25847708-25847708 -v 1"
  decho "$PROG pileup -r $GENOME_DIR/hg19/hg19.fa -i $base/raw_bam/NIC1254A46.bam -q 1 -g chr20:29425717-29425717 -v 1"
}

function biscuittest_somatic_pileup {
  base=MouseWGBS_APCminTumorVsNormal
  [[ -d $base/testresult_somatic_pileup ]] || mkdir $base/testresult_somatic_pileup
  rm -f $base/testresult_somatic_pileup/*
  rname=$(date +%Y_%m_%d)_tumor_vs_normal.vcf
  decho "$PROG pileup -i $base/raw_bam/tumor.bam $base/raw_bam/normal.bam -r $GENOME_DIR/mm10/mm10.fa -o $base/testresult_somatic_pileup/$rname -q 28 -T"
  decho "bgzip $base/testresult_somatic_pileup/rname.gz"
  decho "tabix -p vcf $base/testresult_somatic_pileup/rname.gz"
}

function biscuittest_epiread {
  base=$(pwd)
  [[ -d raw_pileup ]] || cd Smadh3_chr19_chrM
  [[ -d tmp_epiread ]] || mkdir tmp_epiread
  rm -f tmp_epiread/*

  decho "$PROG vcf2bed -t snp raw_pileup/Smadh3_chr19_chrM.vcf.gz >tmp_epiread/Smadh3_chr19_chrM.snp.bed"
  decho "$PROG epiread -n 3 -r $GENOME_DIR/mm10/mm10.fa -i raw_bam/WGBS_Smadh3_chr19_chrM.bam -B tmp_epiread/Smadh3_chr19_chrM.snp.bed -q 20 | gzip -c >tmp_epiread/Smadh3_chr19_chrM.epiread.gz"
  pecho "Collating paired epireads"
  zcat tmp_epiread/Smadh3_chr19_chrM.epiread.gz | sort -k2,2 -k3,3n | awk 'BEGIN{qname="";rec=""}qname==$2{print rec"\t"$5"\t"$6"\t"$7"\t"$8;qname=""}qname!=$2{qname=$2;rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;pair=$3}' | sort -k1,1 -k3,3n | gzip -c >tmp_epiread/Smadh3_chr19_chrM.epiread.paired.gz

  cd $base
}
