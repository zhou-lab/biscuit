************************************
Measure Cytosine Retention and SNP
************************************

**biscuit pileup** computes 1) cytosine retention; and 2) callable SNP mutations.

BISCUIT tries to minimize false positive methylation call and can be stringent in terms of mutations. If there is a read mapped with high confidence and showing a mutation with high base quality score. pileup starts its SNP processing and MAY abolish the methylation calling if the SNP interferes the determination of cytosine retention/conversion.

Feature
#########

- fast, multi-way pileup of multiple samples with VCF output
- computing beta values
- compute genotype, genotype likelihood and genotype quality
- compute somatic score and somatic status when tumor and matched normal samples are provided (-T option)
- de novo infer bisulfite parent strand or from bsmap (ZS) or BWA-meth style tags (YD)
- call ambiguous alternative allele
- distinguish ambiguous alternative allele and multiple alternative allele
- coordinate-sorted output for sorted bam input
- flexible read filtering based on retention number, mapping quality, duplicate and mate pairing.
- flexible base filtering based on base quality, distance to the read ends.
- verbose output of all cytosine context, all the bases, parent strand etc.

Pile up one sample
#####################

To collect methylation and SNP from chr20:1256423-1257433

.. code:: bash

   biscuit pileup -r hg19.fa -i NIC1254A46.bam -q 1 -g chr20:47419734-47419734

::

   chr20   47419734  .   C  .   19   PASS    NS=1;CX=CG;N5=CACGG  
      DP:GT:GP:GQ:SP:CV:BT  12:0/0:1,6,19:19:C7:2:0.00

DP is the total read coverage and CV is the cytosine strand coverage. BT shows beta value (equal to 0.0 in this case). GT, GP and GQ indicate genotype, genotyping likelihood (for 3 genotypes) and genotype quality (for the called genotype). SP shows the allelic support, in this case, all 7 reads after filtering shows C at this particular base. The INFO tag shows the CpG context as well as 5 base sequence environment.

To see how the 12 reads decompose, and why one only gets 0 retention and 2 conversion, one can use `-v 1` to get some extra diagnostic information

::

   DIAGNOSE;
   RN=0;CN=2;
   Bs0=CCCTT; Sta0=66677; Bq0=;7);<; Str0=++---; Pos0=74,68,51,44,23; Rret0=21,19,21,22,18;
   Bs1=CCCCCCC; Sta1=8888888; Bq1=;;):<<2; Str1=+++----; Pos1=69,63,63,62,59,20,1; Rret1=19,19,15,17,19,26,25

The tags are self-explanatory and explained in the VCF header when the verbose tag is given. In brief, the above shows that 2 high quality reads show conversion (CN=2) while no reads show retention (RN=0). 5 reads are of bisulfite Watson (BSW, `ZS:+?` from bsmap, `YD:f` from BWA-meth) strand and 7 reads are on bisulfite Crick (BSC, `ZS:-?` from bsmap, `YD:r` from BWA-meth) strand. Base identities for the BSW reads shows CCCTT (Bs0) and BSC reads shows CCCCCCC (Bs1). Bq0 and Bq1 shows the Phred-scaled base qualities of the bases from these covering reads. Str0 and Str1 shows whether the reads are on forward or reverse strand. Pos0 and Pos1 shows position of the base on the read. Rret shows the number of retention (cytosine methylation) on each read. Sta0 and Sta1 shows the retention-mutation status code internally determined (see below). 
 
Out of the 5 BSW that can inform methylation, 2 reads showing retention is of low base quality and the other is the second last base in the read (a tunable filtering option). On the contrary, 2 other reads suggesting conversion are both high qualities. Hence the output of retention and conversion count can be understood. This example shows how pileup operates on a low coverage, low quality region. One could filter based on the sum of retention and conversion to constrain the beta value computation. Note that the verbose mode `-v 1` also prints positions with no SNP or cytosine methylation. This allows for differentiation of "no mutation" from "no coverage".

code for retention-mutation status
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0: mutation to A;
1: mutation to C;
2: mutation to G;
3: mutation to T;
4: mutation to Y;
5: mutation to R;
6: retention;
7: conversion;
8: reference base;

Pileup multiple samples
###########################

BISCUIT put mutation calling and DNA methylation measurement from multiple samples next to each other when more than one bam is provided as input.

.. code::

   biscuit pileup -i tumor.bam normal.bam -r mm10.fa -g chr19:1-3100000

output

::

   ...
   ##FORMAT=<ID=SC,Number=1,Type=Float,Description="Somatic score">
   ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Variant allele fraction">
   #CHROM  POS  ID  REF   ALT  QUAL    FILTER  INFO    FORMAT  tumor   normal
   chr19 3078956 .  T  A  4   PASS  NS=2  DP:GT:GP:GQ:SP  9:0/0:4,5,33:4:T7,A1  8:0/0:1,7,39:24:T8,A0
   chr19 3078957 .  C  .  21  PASS  NS=2;CX=CHH;N5=TTCTT DP:GT:GP:GQ:SP:CV:BT  8:0/0:1,6,35:21:C7:5:0.00       9:0/0:1,7,39:24:C8:4:0.00
   chr19 3078960 .  C  .  21  PASS  NS=2;CX=CHG;N5=TTCAG DP:GT:GP:GQ:SP:CV:BT  8:0/0:1,6,35:21:C7:5:0.00       10:0/0:1,7,44:27:C9:5:0.00
   ...


Call somatic mutations
^^^^^^^^^^^^^^^^^^^^^^^^^^

If one supply 2 bams to BISCUIT and an option `-T`, BISCUIT would treat the first bam as tumor and the second bam as normal and call somatic mutations.

.. code:: bash

   biscuit pileup -i tumor.bam normal.bam -r mm10.fa -g chr19:1-3100000 -T

output

::

   ...
   chr19 3078956 .  T  A  4   PASS  NS=2;SS=5;SC=0
   DP:GT:GP:GQ:SP  9:0/0:4,5,33:4:T7,A1  8:0/0:1,7,39:24:T8,A0
   ...

outputs somatic states (SS) and somatic score (SC). SS=0,1,2,3,4,5 representing wildtype, germline, somatic, LOH, post-transcriptional modification and unknown respectively; `-x` controls the estimated contamination rate. A higher contamination rate gives more conservative somatic calls (fewer SS=2 calls);

Ambiguous alternative allele
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code::

   pileup -r hg19.fa -i NIC1254A46.bam -q 1 -g chr20:29570686-29570686 -v

outputs

::

   chr20  29570686  .  G  Y  46  PASS  NS=1;CX=CG;N5=ATCGG
      DP:GT:GP:GQ:SP:CV:BT  15:0/1:14,4,37:46:G9,Y4:4:1.00
      DIAGNOSE;RN=4;CN=0;
      Bs0=GGTGTTTGGG;Sta0=8848444888;
      Bq0=68<<=<<<<:;Str0=--++++-++-;Pos0=75,70,55,53,46,41,41,35,22,19;
      Rret0=16,18,21,21,21,23,19,23,21,21;
      Bs1=GGGGG;Sta1=66666;
      Bq1=:*;<<;Str1=++--+;Pos1=62,58,55,34,30;Rret1=27,24,24,22,23

The outputs show alternative allele Y (IUPAC code for `C or T`, supported by 4 reads) when BSC does not suggest alternative allele and there is equal chance of T and C (assuming no prior information of methylation and conversion ratio).

**pileup** differentiates methylation-callable and uncallable (when there is C>T or G>A mutation to confuse methylation calling)

::

   chr20  26138808  .  G  Y  8  PASS  NS=1;CX=CG;N5=TTCGA
      DP:GT:GP:GQ:SP:CV:BT  44:0/1:15,14,143:8:G34,Y6:17:0.76
      DIAGNOSE;RN=13;CN=4;
      Bs0=GGGGGGGGGGGTTGGGGTTTTGGGTG;Sta0=88888888888448888444488848;
      Bq0=:<8<<8<<<<<<:889:9<9988888;Str0=--+++++--++----++--+++-+++;
      Pos0=61,56,53,51,40,38,37,37,27,20,20,20,19,18,14,11,10,10,10,6,6,3,3,2,1,1
      Rret0=4,6,6,6,8,8,8,8,8,10,10,9,9,10,10,11,11,10,10,10,10,12,12,12,11,20;
      Bs1=AAAGGGGGGGGGGGGAAG;Sta1=777666666666666776;
      Bq1=;;;:;;;;<<<<<78;99;Str1=+++--+---++--++---;
      Pos1=75,70,69,67,64,63,62,61,48,45,39,27,19,12,6,5,4,3;
      Rret1=9,10,10,11,12,12,13,13,13,13,13,17,15,16,17,15,15,18

and

::

   chr20  25847708  .  G  Y  20  PASS  NS=1;CX=CG;N5=TCCGT
      DP:GT:GP:GQ:SP:CV:BT  16:0/0:4,8,60:20:G13,Y1:5:1.00
      DIAGNOSE;RN=5;CN=0;
      Bs0=GGGGGGGGGGT;Sta0=88888888884;Bq0=79;;:.<9<8<;Str0=+++-+---+-+;
      Pos0=74,70,65,59,55,50,40,40,36,25,14;Rret0=20,20,23,19,22,13,21,18,22,22,12;
      Bs1=GGGGG;Sta1=66666;Bq1=9:<<:;Str1=--+-+;
      Pos1=69,49,46,37,13;Rret1=18,16,17,17,19

are methylation-callable. But

::

   chr20  29425717  .  G  A  37  PASS  NS=1
      DP:GT:GP:GQ:SP  30:0/0:5,14,111:37:G25,A2
      DIAGNOSE;RN=5;CN=3;
      Bs0=GGTGGGGGGGGGGGGGGAAGGG;Sta0=8848888888888888800888;
      Bq0=8<8;::=;8<<<<<<8*;;:;/;Str0=-+-+-+--+-----+--+++--;
      Pos0=97,76,73,66,65,62,60,59,58,50,47,44,39,37,24,22,17,15,14,8,4,4;
      Rret0=16,17,13,15,15,16,22,17,17,16,15,16,17,17,20,19,19,17,18,19,18,14;
      Bs1=GAGGAAGG;Sta1=67667766;
      Bq1=998;<;9;;Str1=--++-++-;
      Pos1=95,72,61,56,31,27,9,9;Rret1=13,12,11,8,9,10,9,9

is NOT methylation-callable. There is no BT field in FORMAT.
