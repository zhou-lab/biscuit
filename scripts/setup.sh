#!/usr/bin/env bash

# Directory where all assets used in QC.sh live
assets_dir="$1"

# Use <unused> if the file is doesn't exist,
# If <unused>, the corresponding QC section will be skipped
# CpGs
export BISCUIT_CPGBED="${assets_dir}/cpg.bed.gz"
# CpG islands
export BISCUIT_CGIBED="${assets_dir}/cgi.bed.gz"
# Repeat masker bed file
export BISCUIT_RMSK="${assets_dir}/rmsk.bed.gz"
# Merged exon bed file
export BISCUIT_EXON="${assets_dir}/exon.bed.gz"
# Genes
export BISCUIT_GENE="${assets_dir}/genes.bed.gz"
# Locations for the top 100bp bins in GC content
export BISCUIT_TOPGC_BED="${assets_dir}/windows100bp.gc_content.top10p.bed.gz"
# Locations for the bottom 100bp bins in GC content
export BISCUIT_BOTGC_BED="${assets_dir}/windows100bp.gc_content.bot10p.bed.gz"

# QC operations to perform
export BISCUIT_QC_BASECOV=true
export BISCUIT_QC_DUPLICATE=true
export BISCUIT_QC_CPGCOV=true
export BISCUIT_QC_CPGDIST=true
export BISCUIT_QC_UNIFORMITY=true
export BISCUIT_QC_CPGUNIF=true
export BISCUIT_QC_BSCONV=true
export BISCUIT_QC_MAPPING=true
export BISCUIT_QC_BETAS=true
