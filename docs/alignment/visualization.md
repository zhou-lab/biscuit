---
title: "Visualization"
nav_order: 1
parent: Read Mapping
---

The biscuit tview is similar to samtools tview but color reads
considering bisulfite conversion. Make sure you supply the reference
fasta to the command. i.e. 

```bash
$ biscuit tview -g chr19:7525080 input.bam ref.fa
```

![tview screenshot](/biscuit/assets/2017_05_02_biscuit_tview_figure.png)

For reference, CpGs are marked in red, C and Gs in other context are
marked in blue. For reads, retended C/Gs are marked in red, converted
C/Gs are marked in blue. Mismatches not related to bisulfite
conversion are marked in yellow.
