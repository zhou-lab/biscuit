# BISCUIT [![Travis-CI Build Status](https://travis-ci.org/huishenlab/biscuit.svg?branch=master)](https://travis-ci.org/huishenlab/biscuit) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.48262.svg)](http://dx.doi.org/10.5281/zenodo.48262)

BISulfite-seq CUI Toolkit (BISCUIT) is a utility for analyzing sodium bisulfite
conversion-based DNA methylation/modification data. It was written to perform
alignment, DNA methylation and mutation calling, and allele specific methylation
from bisulfite sequencing data.

# Download and Install

Instructions for downloading and installing BISCUIT can be found in the User
Guide [Download and Install](https://huishenlab.github.io/biscuit/#download-and-install)
page.

All releases of BISCUIT are available on
[GitHub](https://github.com/huishenlab/biscuit/releases).

Note, to use the BISCUIT QC script, the following tools must be installed and in
the user's `PATH` environment variable:

  - BISCUIT
  - [samtools](http://www.htslib.org/)
  - [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
  - [GNU awk](https://www.gnu.org/software/gawk/manual/gawk.html)
  - [GNU parallel](https://www.gnu.org/software/parallel/)

It is also useful to have [samblaster](https://github.com/GregoryFaust/samblaster)
installed for marking duplicate reads during the alignment phase.

# Usage and Documentation

A User Guide has been created to provide useful documentation regarding the
usage of BISCUIT. The guide can be found at:
[https://huishenlab.github.io/biscuit/](https://huishenlab.github.io/biscuit/).

# Issues

Any issues with BISCUIT can be submitted on the Issues page:
[https://github.com/huishenlab/biscuit/issues](https://github.com/huishenlab/biscuit/issues).

<!-- # Links -->

<!--   - BISCUIT Publication -->
<!-- > Link to publication -->

# Acknowledgements

 * lib/aln was adapted from Heng Li's BWA-mem code.
 * lib/htslib was subtree-ed from the htslib library.
 * lib/klib was subtree-ed from Heng Li's klib.
