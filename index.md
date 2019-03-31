---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults
title: HOME
layout: default
description: "BISulfite-seq CUI Toolkit (BISCUIT) BISulfite-seq CUI Toolkit
  (BISCUIT) is a utility suite for analyzing sodium bisulfite
  conversion-based DNA methylation/modification data. It was written
  to perform alignment, DNA methylation and mutation calling, allele
  specific methylation from bisulfite sequencing data."
nav_order: 1
---

# BISCUIT - Focus on Sequencing Data with Bisulfite Conversion
{: .fs-9 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/zwdzwd/biscuit){: .btn .fs-5 .mb-4 .mb-md-0 }

## Getting Started

To start, just download the precompiled versions from [Github](https://github.com/zwdzwd/biscuit/releases/latest). 
One can do this in terminal using the following one-liner:

For mac OS,
```bash
$ curl -OL $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest |
    grep browser_download_url | grep macOS | cut -d '"' -f 4) && chmod a+x biscuit*
```

For linux,
```bash
$ curl -OL $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest | 
    grep browser_download_url | grep x86_64 | cut -d '"' -f 4) && chmod a+x biscuit*
```

## Source code
All releases are available [here](https://github.com/zwdzwd/biscuit/releases/). Note after v0.2.0, if you choose to use Git, make sure use `git clone --recursive` to get the submodules.

You can also download with Curl using
```bash
$ curl -O $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest | grep browser_download_url | grep release-source.zip | cut -d '"' -f 4)
```

## Compile

If you download source code, BISCUIT can be compiled through,

```bash
$ unzip release.zip
$ cd biscuit-release
$ make
```

The created `biscuit` binary is the main entry point.

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
$ biscuit align -t 10 GRCh38.fa fastq1.fq.gz fastq2.fq.gz | samtools sort -T . -O bam -o output.bam
$ samtools index output.bam
$ samtools flagstat output.bam >output.bam.flagstat
```

See [here](https://github.com/zwdzwd/biscuit/wiki/Measure-cytosine-retention-and-SNP) for more information.


## Visualize alignment

The `tview` subroutine colors the alignments in bisulfite mode. [Here](https://github.com/zwdzwd/biscuit/wiki/Visualize-reads-with-bisulfite-conversion) is a screenshot.

```bash
$ biscuit tview -g chr19:7525080 input.bam ref.fa
```
Unlike samtools, in this subroutine, a reference fasta file is mandatory so that bisulfite conversion can be identified.

## Mark duplicate reads

This step is optional. The mark duplicate of BISCUIT is bisulfite strand aware.
```bash
$ biscuit markdup input.bam output.bam
```

## Pileup

Like samtools, BISCUIT extract DNA methylation as well as genetic information. The following shows how to produce a tabix-indexed VCF file.
```bash
$ biscuit pileup GRCh38.fa input.bam -o output.vcf -q 20
$ bgzip output.vcf
$ tabix -p vcf output.vcf.gz
```

## Make bed files

The following extract CpG beta values from the VCF file.
```bash
$ biscuit vcf2bed -k 10 -t cg input.vcf.gz
```

`-t` can also take

  * `snp` - SNP information
  * `c` - all cytosines
  * `hcg` - HCG for NOMe-seq
  * `gch` - GCH for NOMe-seq
  
`-e` output sequence contexts.

## QC your run

A bash script is provided to simplify QC procedure.
```bash
$ ./scripts/QC.sh -v input.vcf setup_file sample.name input.bam
```
generates QC information that can be picked up by MultiQC. Setup files for [hg19](http://zwdzwd.io/BISCUITqc/hg19_QC_assets.zip), [hg38](http://zwdzwd.io/BISCUITqc/hg38_QC_assets.zip) and [mm10](http://zwdzwd.io/BISCUITqc/mm10_QC_assets.zip) are provided. Other genome builds can be built following the same format.

## EPI-reads and allele-specific methylation

Following illustrates how to produce `epiread` which carries the information of epi-haplotype.
```bash
$ biscuit epiread -r GRCh38.fa -i input.bam -B snp.bed
```

To test all SNP-CpG pair,
```bash
$ biscuit epiread -r GRCh38.fa -P -i input.bam -B snp.bed
```
Details can be found [here](https://github.com/zwdzwd/biscuit/wiki/Convert-to-epiread-format).

```bash
sort -k1,1 -k2,2n -k3,3n in.epiread >out.epiread
biscuit asm out.epiread >out.asm
```

## Validate bisulfite conversion label

Sometimes, the bisulfite conversion labels in a given alignment are inaccurate, conflicting or ambiguous. The `bsstrand` command summarizes these labels given the number of C>T, G>A substitutions. It can correct inaccurate labels as an option.
```bash
$ biscuit bsstrand -g "chr1:1000000-1050000" GRCh37.fa input.bam 
```
gives something like
```
Mapped reads: 16688
Unmapped reads: 29
Corrected reads: 0 (0.00%)

Confusion counts:
orig\infer     BSW (f)      BSC (r)      Conflict (c) Unknown (u)
     BSW (f):   8426         42           4            12
     BSC (r):   15           8167         3            19
Conflict (c):   0            0            0            0
 Unknown (u):   0            0            0            0
```


The inferred `YD` tag gives the following
- f: foward/Waston strand
- r: reverse/Crick strand
- c: conflicting strand information
- u: unintelligible strand source (unknown)

`YD` is inferred based on the count of `C>T` (`nCT`) and `G>A` (`nGA`) observations in each read according to the following rule: If both `nCT` and `nGA` are zero, `YD = u`. `s = min(nG2A,nC2T)/max(nG2A,nC2T)`. if `nC2T > nG2A && (nG2A == 0 || s <= 0.5)`, then `YD = f`. if `nC2T < nG2A && (nC2T == 0 || s <= 0.5)`, then `YD = r`. All other scenarios gives `YD = c`. `-y` append `nCT`(YC tag) and `nGA`(YG tag) in the output bam.

## Summarize and filter reads by bisulfite conversion

```bash
$ biscuit bsconv -g "chr1:1000000-1050000" GRCh37.fa input.bam
```
For some library preparation, incomplete conversion are enriched in a subset of reads that needs to be filtered. This command transforms bam into one that contains file `ZN` tag e.g., `ZN:Z:CA_R0C11,CC_R1C14,CG_R0C2,CT_R1C5`. This tag summarizes counts of retention and conversion for four different cytosine contexts `CpA`, `CpC`, `CpG` and `CpT`. It contains a minimum threshold of `CpA`, `CpC`, `CpT` or `CpH` in general. The `-b` option outputs the summary in tables instead of as tags in the BAM file.


## About the project

This package is made by the folks from Van Andel Research Institute
with help from prior code base from the internet.

### License

### Contributing
