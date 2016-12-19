# BISCUIT [![Travis-CI Build Status](https://travis-ci.org/zwdzwd/biscuit.svg?branch=master)](https://travis-ci.org/zwdzwd/biscuit) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.48262.svg)](http://dx.doi.org/10.5281/zenodo.48262)

BISulfite CUI Toolkit (BISCUIT) is a utility suite for analyzing sodium bisulfite conversion-based DNA methylation/modification data. It was written to perform alignment, DNA methylation and mutation calling, allele specific methylation from bisulfite sequencing data.

# Download and Install

To install BISCUIT,

```bash
$ wget https://github.com/zwdzwd/biscuit/archive/master.zip
$ unzip master.zip
$ cd biscuit-master
$ make
```

The created `biscuit` binary is the main entry point.

Stable releases can be downloaded [here](https://github.com/zwdzwd/biscuit/releases/).

<!-- User Guide is available [here](https://github.com/zwdzwd/biscuit/wiki). -->

# Usage

## Index reference for alignment

```bash
biscuit index GRCh38.fa
```
The index of BISCUIT composed of the 2-bit packed reference (`.bis.pac`, `.bis.amb`, `.bis.ann`). The suffix array and FM-index of the parent strand (`.par.bwt` and `.par.sa`) and the daughter strand (`.dau.bwt` and `.dau.sa`).

## Read alignment

The following snippet shows how BISCUIT can be used in conjunction with [samtools](https://github.com/samtools/samtools) to produce indexed alignment BAM file. 
```bash
$ biscuit align GRCh38.fa -t 10 fastq1.fq.gz fastq2.fq.gz | samtools sort -T . -O bam -o output.bam
$ samtools index output.bam
$ samtools flagstat output.bam >output.bam.flagstat
```

## Mark duplicate reads

This step is optional. The mark duplicate of BISCUIT is bisulfite strand aware.
```bash
$ biscuit markdup input.bam output.bam
```

## Pileup

Like samtools, BISCUIT extract DNA methylation as well as genetic information. The following shows how to produce a tabix-indexed VCF file.
```bash
$ biscuit pileup -r GRCh38.fa -i input.bam -o output.vcf -q 20
$ bgzip output.vcf
$ tabix -p vcf output.vcf.gz
```

## Make bed files

The following extract CpG beta values from the VCF file.
```bash
$ biscuit vcf2bed -k 10 -t cg input.vcf.gz
```

# Acknowledgements

 * lib/aln was adapted from Heng Li's BWA-mem code.
 * lib/htslib was subtree-ed from the htslib library.
 * lib/klib was subtree-ed from Heng Li's klib.
