# BISCUIT

BISulfite-seq CUI Toolkit (BISCUIT) is a utility for analyzing sodium bisulfite
conversion-based DNA methylation/modification data. It was written to perform
alignment, DNA methylation and mutation calling, and allele specific methylation
from bisulfite sequencing data.

## Publication

If you use BISCUIT, kindly [cite](https://doi.org/10.1093/nar/gkae097):

```
Wanding Zhou, Benjamin K Johnson, Jacob Morrison, Ian Beddows,
James Eapen, Efrat Katsman, Ayush Semwal, Walid Abi Habib, Lyong Heo,
Peter W Laird, Benjamin P Berman, Timothy J Triche, Hui Shen,
BISCUIT: an efficient, standards-compliant tool suite for simultaneous
    genetic and epigenetic inference in bulk and single-cell studies,
Nucleic Acids Research, Volume 52, Issue 6, 12 April 2024, Page e32,
https://doi.org/10.1093/nar/gkae097
```

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

It is also useful to have [dupsifter](https://github.com/huishenlab/dupsifter)
installed for marking duplicate reads during the alignment phase.

# Usage and Documentation

A User Guide has been created to provide useful documentation regarding the
usage of BISCUIT. The guide can be found at:
[https://huishenlab.github.io/biscuit/](https://huishenlab.github.io/biscuit/).

# Issues

Any issues with BISCUIT can be submitted on the Issues page:
[https://github.com/huishenlab/biscuit/issues](https://github.com/huishenlab/biscuit/issues).

# Acknowledgements

 * lib/aln was adapted from Heng Li's BWA-mem code.
 * This work is supported by NIH/NCI R37CA230748 and U24CA210969.
