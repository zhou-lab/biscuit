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

# BISCUIT - Decipher Sequencing Data with Bisulfite Conversion
{: .fs-9 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/zwdzwd/biscuit){: .btn .fs-5 .mb-4 .mb-md-0 }

### Getting Started

To start, just download the [Precompiled Binaries](https://github.com/zwdzwd/biscuit/releases/latest).
from Github (currently only supports latest versions of Linux and MacOSX).
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

### Compile from Source Code

You can compile from source code using the following command. Note if
you choose to clone from Github, make sure you specify `git clone --recursive`
to get the submodules.

```bash
$ curl -OL $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest | 
    grep browser_download_url | grep release-source.zip | cut -d '"' -f 4)
$ unzip release.zip
$ cd biscuit-release
$ make
```

The created `biscuit` binary is the main entry point.

### Usage

Functionalities are listed by just typing `biscuit` in terminal.

```bash
$ biscuit
```
```
Program: BISCUIT (BISulfite-seq CUI Toolkit)
Version: 0.3.9.20190330
Contact: Wanding Zhou <wanding.zhou@vai.org>

Usage:   biscuit <command> [options]

Command:
  -- Read mapping
     index         index reference genome sequences in the FASTA format
     align         align bisulfite treated short reads using adapted BWA-mem algorithm

  -- BAM operation
     tview         text alignment viewer with bisulfite coloring
     markdup       mark duplicates on the same bisulfite strand
     bsstrand      validate/correct bisulfite conversion strand label (YD tag)
     bsconv        summarize/filter reads by bisulfite conversion (ZN tag)
     cinread       print cytosine-read pair in a long form.

  -- Base summary
     pileup        pileup cytosine and mutations.
     vcf2bed       convert VCF to bed graph.
     mergecg       merge C and G in CpG context.

  -- Epireads
     epiread       convert bam to epiread format
     rectangle     convert epiread to rectangle format
     asm           test allele specific methylation
```

### About the project

This package is made by the folks from Van Andel Research Institute
with help from prior code base from the internet.

### Acknowledgement

 - lib/aln was adapted from Heng Li's BWA-mem code.
 - lib/htslib was submoduled from the htslib library.
 - lib/klib was submoduled from Heng Li's klib.

### Reference

In preparation
