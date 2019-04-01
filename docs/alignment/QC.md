---
title: Quality Control
nav_order: 2
parent: Read Mapping
---
# Quality Control

A bash script is provided to simplify QC procedure.
```bash
$ ./scripts/QC.sh -v input.vcf setup_file sample.name input.bam
```
generates QC information that can be picked up by MultiQC. Setup files for [hg19](http://zwdzwd.io/BISCUITqc/hg19_QC_assets.zip), [hg38](http://zwdzwd.io/BISCUITqc/hg38_QC_assets.zip) and [mm10](http://zwdzwd.io/BISCUITqc/mm10_QC_assets.zip) are provided. Other genome builds can be built following the same format.

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
