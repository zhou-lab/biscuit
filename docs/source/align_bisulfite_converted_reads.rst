*************************************
Align bisulfite-converted reads
*************************************

biscuit index and biscuit align adapt the bwa index and bwa mem for bisulfite-treated short reads

Design
#########

- assymmetric scoring for C to T and G to A in mapping.
- produce consistent mapping quality calculation with destination strand specification.
- produce consistent NM and MD tags under assymmetric scoring.
- produce ZR and ZC tags for retention count and conversion count
- separate seeding for parent and daughter strands for mapping efficiency
- disk space economic indices with no need to store a bisulfite converted reference.
- BWA-mem parameters visible to the users.
- no separate installation of BWA. Dependencies easily met.
- robust to OS build, read processing and reference processing occur within computing threads of a single process.
- support for both single- and paired-end reads
- Optional parent strand and daughter strand restriction for both single- and paired-end reads.
- Optional BSW/top/BSC/bottom strand restriction, tightly integrated in mapping.

Force mapping to a particular strand
#####################################

Force mapping to the parent strand
-B 3

Force mapping to the daughter strand (mapping probes, for example)
-B  1
