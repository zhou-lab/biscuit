---
title: Extract Methylation
nav_order: 4
---

# Extract Methylation



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

### Merge neighboring C and G in CpG context
