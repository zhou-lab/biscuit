
# bin/biscuit align ~/reference/mm10/mm10.fa /data/largeS2/pl-bs/2015-06-09-mouse-WGBS/run3_H7CVYBGXX/data/fastq/PL4-10-15WGBS1_L000_R1_001.fastq.gz /data/largeS2/pl-bs/2015-06-09-mouse-WGBS/run3_H7CVYBGXX/data/fastq/PL4-10-15WGBS1_L000_R2_001.fastq.gz -t 70 | samtools view -bS - > /data/largeS2/pl-bs/H7CVYBGXX_WGBS1.bam

PROG=~/tools/biscuit/development/biscuit/bin/biscuit
GENOME_DIR=~/references/

function pecho() {
  echo "[$(date)] Running:  "$@ >&2
}

function decho() {
  echo "[$(date)] Running:  "$@ >&2
  eval $@
}

function biscuittest_pileup {
  base=$(pwd)
  [[ -d raw_bam ]] || cd Smadh3_chr19_chrM
  [[ -d tmp_pileup ]] || mkdir tmp_pileup
  rm -f tmp_pileup/*
  decho "$PROG pileup -i raw_bam/WGBS_Smadh3_chr19_chrM.bam -n 3 -r $GENOME_DIR/mm10/mm10.fa -o tmp_pileup/Smadh3_chr19_chrM.vcf -q 28"
  decho "bgzip tmp_pileup/Smadh3_chr19_chrM.vcf"
  decho "tabix -p vcf tmp_pileup/Smadh3_chr19_chrM.vcf.gz"
  cd $base
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
