---
layout: default
title: Working with Epialleles
permalink: /epiread/
---

## Generating epireads

Epiread format is a compact way of storing CpG retention pattern as well as SNP information on the same read. It is useful for estimating epiallele fraction and clonal structure/cell population (Li et al. Genome Biology 2014, Zheng et al. Genome Biology 2014). BISCUIT extends the original epiread format proposed by the methpipe team by allowing SNP information displayed. The columns indicate: **1)** chromosome name; **2)** read name; **3)** read position in paired-end sequencing; **4)** bisulfite strand (bisulfite Watson or Crick); **5)** position of the cytosine in the first CpG (0-based); **6)** retention pattern ("C" for retention and "T" for conversion) for all CpGs covered; **7)** position of the first SNP if SNP location file is provided; **8)** base call of all SNPs covered.
```
chr19   NS500653:8:HF5FGBGXX:3:12402:11299:9856 1       +       3040315 CCCCTCCC        .       .
chr19   NS500653:8:HF5FGBGXX:3:12609:17196:5738 1       +       3055491 T       .       .
chr19   NS500653:8:HF5FGBGXX:1:23304:20253:14257        1       +       3078472 CC      3078510 T
chr19   NS500653:8:HF5FGBGXX:4:22509:19067:5776 1       +       3078472 CC      3078510 C
chr19   NS500653:8:HF5FGBGXX:4:22611:9688:19912 2       +       3078946 C       .       .
chr19   NS500653:8:HF5FGBGXX:4:22602:25913:14920        2       +       3078946 CC      3078982 A
```

To produce an epiread format, one needs to run the `epiread` subcommand with a bam file, a fasta file for the reference and optionally a bed file for SNPs.
```bash
biscuit epiread -r genome.fa -i alignment.bam [-B snps.bed]
```
The SNP bed file can be obtained through `biscuit vcf2bed -t snp pileup.vcf.gz`. If no SNP file is supplied, the output does not include extra columns for SNP.
`cut -f 1,5,6` gets the original epiread format.

## Paired-end epireads

DNA methylation information from both mate reads are physically "phased" molecular events that can be traced back to the same DNA molecule. The read name and read position in the single-end epiread output can be used to collate mate reads in a read pair. The default behavior of epiread subcommand focuses only on the primary mapping. The following awk gives a nice and compact file for paired-end epiread format:
```bash
sort -k2,2 -k3,3n singleend_epiread |
  awk 'BEGIN{qname="";rec=""}
       qname==$2{print rec"\t"$5"\t"$6"\t"$7"\t"$8;qname=""}
       qname!=$2{qname=$2;rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;pair=$3}'
```
In the following output, columns represent: **1)** chromosome; **2)** bisulfite strand; **3)** location of the cytosine in the first CpG covered (0-based) in read 1; **4)** retention pattern of all the CpG covered in read 1; **5)** location of the first SNP covered in read 1; **6)** base call of all the SNPs covered in read 1; **7)** location of the cytosine in the first CpG covered (0-based) in read 2; **8)** retention pattern of all the CpG covered in read 2; **9)** location of the first SNP covered in read 2; **10)** base call of all the SNPs covered in read 2;
```
chr19   -       3083513 CCCCCCC 3083495 ATT     3083513 CCCCCCC 3083495 ATT
chr19   -       3083545 CCTCCCCCCT      .       .       3083527 CCCCTCCCCCC     3083523 AG
chr19   +       3083616 TTTTTTTTTTT     .       .       3083722 TTTTTTTTTTT     .       .
chr19   -       3083630 CCCCCTCCCCCCT   .       .       3083616 TCCCCCTCCCC     .       .
chr19   +       3083638 TTTTTTTTTTTT    .       .       3083705 TTTTTTTTTTT     .       .
```

## NOMe-seq epireads

The `-N` option allows listing GCH and HCG retention states in NOMe-seq data be listed side by side. E.g.,
<!--
$$$ mkdir -p test/NOMeSeq_HCT116_chr18_chr19_chrM/out_epiread
@@ biscuit : biscuit-develop
@@ HCT116.vcf.gz : test/NOMeSeq_HCT116_chr18_chr19_chrM/raw_pileup/HCT116_chr18_chr19_chrM.vcf.gz
@@ hg19.fa : /home/wanding.zhou/references/hg19/hg19.fa
@@ HCT116.bam : test/NOMeSeq_HCT116_chr18_chr19_chrM/raw_bam/HCT116_chr18_chr19_chrM.bam
@@ snp.bed : test/NOMeSeq_HCT116_chr18_chr19_chrM/out_epiread/snp.bed
@@ epiread.gz : test/NOMeSeq_HCT116_chr18_chr19_chrM/out_epiread/epiread.gz
@@ epiread_paired.gz : test/NOMeSeq_HCT116_chr18_chr19_chrM/out_epiread/epiread_paired.gz
-->
```bash
$$+ biscuit vcf2bed -t snp HCT116.vcf.gz >snp.bed
$$+ biscuit epiread -n 3 -N -r hg19.fa -i HCT116.bam -B snp.bed -N -q 20 | gzip -c >epiread.gz
# collating paired epireads
$$+ zcat epiread.gz | sort -k2,2 -k3,3n | awk 'BEGIN{qname="";rec=""}qname==$2{print rec"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10;qname=""}qname!=$2{qname=$2;rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10;pair=$3}' | sort -k1,1 -k3,3n | gzip -c >epiread_paired.gz
```
<!--
##compare test/NOMeSeq_HCT116_chr18_chr19_chrM/out_epiread/snp.bed vs test/NOMeSeq_HCT116_chr18_chr19_chrM/golden_epiread/snp.bed
##compare test/NOMeSeq_HCT116_chr18_chr19_chrM/out_epiread/epiread.gz vs test/NOMeSeq_HCT116_chr18_chr19_chrM/golden_epiread/epiread.gz
##compare test/NOMeSeq_HCT116_chr18_chr19_chrM/out_epiread/epiread_paired.gz vs test/NOMeSeq_HCT116_chr18_chr19_chrM/golden_epiread/epiread_paired.gz
$+ epiread_paired.gz
-->
the paired output looks like,
```text
chr18  +  10689  CCC  10663  TTCT  10696  T  10689  CCC  10663  TTCT  10696  T
chr18  +  10689  CCC  10694  TT  10696  T  10689  CCC  10694  TT  10696  T
chr18  -  10689  CCC  10699  TTT  10696  T  10689  CCC  10699  TTTT  10696  T
chr18  +  10703  CCCCCC  10694  CCCCCTCTTTT  10696  T  10753  CCCCC  10743  CTCTTTTT
    .  .
chr18  -  11134  CCCC  11131  CTTT  .  .  11134  CCCC  11131  CTTT  .  .
chr18  -  11344  CC  11339  CCTTT  .  .  11344  CC  11339  CCTT  .  .
```
The columns are: **1)** chromosome; **2)** bisulfite strand; **3)** location of the cytosine in the first CpG covered (0-based) in read 1; **4)** retention pattern of all the CpG covered in read 1; **5)** location of the cytosine in the first GpCpH covered (0-based) in read 1; **6)** retention pattern of all the GpCpH covered in read 1; **7)** location of the first SNP covered in read 1; **8)** base call of all the SNPs covered in read 1; **9)** location of the cytosine in the first CpG covered (0-based) in read 2; **10)** retention pattern of all the CpG covered in read 2; **11)** location of the cytosine in the first GpCpH covered (0-based) in read 2; **12)** retention pattern of all the GpCpH covered in read 2; **13)** location of the first SNP covered in read 2; **14)** base call of all the SNPs covered in read 2; To be exact, BISCUIT only counts HpCpG for CpG. But BISCUIT records the location of C regardless whether it is the G or the C that is measured. In other words, if there is a TGCGA and a read was bisulfite-treated on the Crick strand, BISCUIT only records the 4th G but not the 3rd C but still use the position of the 3rd C as the position of the CpG. This is to be consistent with standard BS-seq.

