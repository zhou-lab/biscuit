# BISCUIT [![Travis-CI Build Status](https://travis-ci.org/zwdzwd/biscuit.svg?branch=master)](https://travis-ci.org/zwdzwd/biscuit) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.48262.svg)](http://dx.doi.org/10.5281/zenodo.48262)
---

BISulfite CUI Toolkit (BISCUIT) is a utility suite for analyzing sodium bisulfite conversion-based DNA methylation/modification data. It was written to perform alignment, DNA methylation and mutation calling, allele specific methylation from bisulfite sequencing data.

# Download and Install

The latest release can be downloaded [here](https://github.com/zwdzwd/biscuit/releases/latest)

To install BISCUIT,

```bash
$ wget https://github.com/zwdzwd/biscuit/archive/master.zip
$ unzip biscuit-master.zip
$ cd biscuit-master
$ make
```

This create BISCUIT binary.


User Guide is available [here](https://github.com/zwdzwd/biscuit/wiki).

# Acknowledgements

 * lib/aln was adapted from Heng Li's BWA-mem code.
 * lib/htslib was subtree-ed from the htslib library.
 * lib/klib was subtree-ed from Heng Li's klib.
