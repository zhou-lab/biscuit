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

# BISCUIT - Understand Sequencing Data with Bisulfite Conversion
{: .fs-9 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/zwdzwd/biscuit){: .btn .fs-5 .mb-4 .mb-md-0 }

---

### Getting Started

To start, just download the [Precompiled Binaries](https://github.com/zwdzwd/biscuit/releases/latest).
from Github (currently only supports latest versions of Linux and MacOSX).
One can do this in terminal using the following one-liner:

For mac OS,
```bash
$ curl -o biscuit \
    -L $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest |
    grep browser_download_url | grep darwin | cut -d '"' -f 4) && chmod a+x biscuit
```

For linux,
```bash
$ curl -o biscuit \
    -OL $(curl -s https://api.github.com/repos/zwdzwd/biscuit/releases/latest | 
    grep browser_download_url | grep linux | cut -d '"' -f 4) && chmod a+x biscuit
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

### Overview of Functionalities

See the following list for an overview of different
functionalities/subcommands provided in `biscuit`. You can also see
this by just typing `biscuit` in terminal.

#### Read Mapping

  - `index` Index reference genome (see 
    [Read Mapping]({{ site.baseurl }}{% link docs/alignment/alignment.md %}))
  - `align` Map bisulfite converted short reads to reference (see
    [Read Mapping]({{ site.baseurl }}{% link docs/alignment/alignment.md %}))

#### BAM Operation

  - `tview` View read mapping in terminal with bisulfite coloring (see
    [Visualization]({{ site.baseurl }}
    {% link docs/alignment/visualization.md %}))
  - `markdup` Mark read duplication (see [BAM operation]({{ site.baseurl }}
    {% link docs/alignment/bam_operation.md %}))
  - `bsstrand` Investigate bisulfite conversion strand label (see
    [Quality Control]({{ site.baseurl }}
    {% link docs/alignment/QC.md %}))
  - `bsconv` Investigate bisulfite conversion rate (see
    [Quality Control]({{ site.baseurl }}
    {% link docs/alignment/QC.md %}))
  - `cinread` Print cytosine-read pair in a long form (see
    [BAM operation]({ site.baseurl }}
    {% link docs/alignment/bam_operation.md %}))

#### Methylation, SNP Extraction

  - `pileup` Generate standard-compliant VCF (see 
    [Read Pileup]({{ site.baseurl }}
    {% link docs/pileup.md %}))
  - `vcf2bed` Extract mutation, methylation from VCF.
    (see [Extract Methylation]({{ site.baseurl }}
    {% link docs/methylextraction.md %}) and 
    [Extract Mutation]({{ site.baseurl }}
    {% link docs/methylextraction.md %}))
  - `mergecg` Merge neighboring C and G in CpG context.
    (see [Extract Methylation]({{ site.baseurl }}
    {% link docs/methylextraction.md %}))
  
#### Epi-read & Epi-allele

  - `epiread` Convert bam to epi-read format (see
    [Epi-read & Epi-allele]({{ site.baseurl }}
    {% link docs/Epiread.md %}))
  - `rectangle` Convert epi-read to rectangle format (see
    [Epi-read & Epi-allele]({{ site.baseurl }}
    {% link docs/Epiread.md %}))
  - `asm` Test allele-specific methylation. (see
    [Allele-specific Methylation]({{ site.baseurl }}
    {% link docs/allele_meth.md %}))

### About the project

This package is made by the folks from Van Andel Research Institute
with help from prior code base from the internet.

### Acknowledgement

 - lib/aln was adapted from Heng Li's BWA-mem code.
 - lib/htslib was submoduled from the htslib library.
 - lib/klib was submoduled from Heng Li's klib.

### Reference

In preparation
