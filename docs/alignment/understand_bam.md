---
title: Understand BAM
parent: Read Mapping
nav_order: 1
---
# Understand BISCUIT bams

The BISCUIT bam is slightly different from usual BAM files in tags and
how insert size are computed.

### TAGS

- `NM` - number of mismatches. This excludes cytosine conversions.
- `XA` - location of the suboptimal alignments
- `XB` - integer pair. The first in pair indicates the number of
  suboptimal mapping in the primary/non-decoy chromosomes. The second
  in pair indicates the number of suboptimal mapping in the ALT/decoy
  chromosomes. E.g., `10,5` means 10 suboptimal alignment exists on
  primary/non-decoy chromosomes and 5 exists on ALT/decoy chromosomes.
- `ZC` - number of cytosine conversion
- `ZR` - number of cytosine retention
- `AS` - best alignment score
- `XS` - suboptimal alignment score. This is usually equal or under
  AS. In rare cases, pairing would cause a `XS` higher than `AS`
- `MD` - location of the mismatches, following samtools convention.
- `PA` - ratio of score (AS) / alt_score (XS), the higher the ratio,
  the more accurate the position
- `SA` - other parts of a chimeric primary mapping
- `YD` - Bisulfite conversion strand label, `f` for forward/Watson and
  `r` for reverse/Crick, a la BWA-meth.

See also [BAM Operation]({{ site.baseurl }}{% link docs/alignment/bam_operation.md %})
for how to add `ZN` tag for cytosine conversion under `CpG`, and other
non-CpG dinucleotide sequence contexts.

### Insert Size

Insert size/TLEN is printed in different ways from BWA. Here the tlen
is the actual insert size:

```
reverse mate read's right-most coordinate - forward mate read's left-most coordinate.
```


