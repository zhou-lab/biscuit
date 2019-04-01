---
title: Read Pileup
nav_order: 3
---

# Generate Standard VCF Output

### The pileup subcommand

Like samtools, BISCUIT extract DNA methylation as well as genetic information. The following shows how to produce a tabix-indexed VCF file.
```bash
$ biscuit pileup GRCh38.fa input.bam -o output.vcf -q 20
$ bgzip output.vcf
$ tabix -p vcf output.vcf.gz
```
