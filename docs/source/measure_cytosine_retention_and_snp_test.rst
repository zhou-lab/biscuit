*****************************************
pileup human TCGA WGBS NIC1254A46 focal
*****************************************

<!---
$$$ mkdir -p test/NIC1254A46/out_pileup
@@ biscuit : biscuit-develop
@@ NIC1254A46.bam : test/NIC1254A46/raw_bam/NIC1254A46.bam
@@ hg19.fa : /home/wanding.zhou/references/hg19/hg19.fa
@@ out.vcf : test/NIC1254A46/out_pileup/pileup_chr20_29570686.vcf
-->
```bash
$ biscuit pileup -i NIC1254A46.bam -r hg19.fa -g chr20:29570686-29570686
```
returns
```text
chr20  29570686  .  G  Y  46  PASS  NS=1;CX=CG;N5=ATCGG  DP:GT:GP:GQ:SP:CV:BT  15:0/1:14,4,37:46:G9,Y4:4:1.00
```

```bash
$ biscuit pileup -i NIC1254A46.bam -r hg19.fa -g chr20:26138808-26138808 -v 1 -n 999
```
returns
```text
chr20  26138808  .  G  Y  5  PASS  NS=1;CX=CG;N5=TTCGA  DP:GT:GP:GQ:SP:CV:BT  45:0/1:15,14,147:5:G35,Y6:18:0.72
    DIAGNOSE;RN=13;CN=5;Bs0=GGGGGGGGGGGTTGGGGTTTTGGGTG;Sta0=88888888888448888444488848;Bq0=:<8<<8<<<<<<:889:9<9988888;Str0=--+++++--++----++--+++-+++;Pos0=61,56,53,51,40,38,37,37,27,20,20
```

```bash
$ biscuit pileup -i NIC1254A46.bam -r hg19.fa -g chr20:25847708-25847708 -v 1
```
returns
```text
chr20  25847708  .  G  Y  20  PASS  NS=1;CX=CG;N5=TCCGT  DP:GT:GP:GQ:SP:CV:BT  16:0/0:4,8,60:20:G13,Y1:5:1.00
    DIAGNOSE;RN=5;CN=0;Bs0=GGGGGGGGGGT;Sta0=88888888884;Bq0=79;;:.<9<8<;Str0=+++-+---+-+;Pos0=74,70,65,59,55,50,40,40,36,25,14;Rret0=20,20,23,19,22,13,21,18,22,22,12;Bs1=GGGGG;Sta1=66666;
```

```bash
$ biscuit pileup -i NIC1254A46.bam -r hg19.fa -g chr20:29425717-29425717 -v 1
```
returns
```text
chr20  29425717  .  G  A  37  PASS  NS=1  DP:GT:GP:GQ:SP  30:0/0:5,14,111:37:G25,A2
    DIAGNOSE;RN=5;CN=3;Bs0=GGTGGGGGGGGGGGGGGAAGGG;Sta0=8848888888888888800888;Bq0=8<8;::=;8<<<<<<8*;;:;/;Str0=-+-+-+--+-----+--+++--;Pos0=97,76,73,66,65,62,60,59,58,50,47,44,39,37,24,22,1
```

# somatic pileup

<!---
$$$ mkdir -p test/MouseWGBS_APCminTumorVsNormal/out_pileup
@@ biscuit : biscuit-develop
@@ tumor.bam : test/MouseWGBS_APCminTumorVsNormal/raw_bam/tumor.bam
@@ normal.bam : test/MouseWGBS_APCminTumorVsNormal/raw_bam/normal.bam
@@ mm10.fa : /home/wanding.zhou/references/mm10/mm10.fa
@@ mouseAPC.vcf : test/MouseWGBS_APCminTumorVsNormal/out_pileup/pileup.vcf
@@ mouseAPC.vcf.gz : test/MouseWGBS_APCminTumorVsNormal/out_pileup/pileup.vcf.gz
-->
```bash
$$ biscuit pileup -i tumor.bam normal.bam -r mm10.fa -g chr19:1-3100000 -o mouseAPC.vcf -T
$$ bgzip -f mouseAPC.vcf
$$ tabix -p vcf mouseAPC.vcf.gz
```
<!---
##compare test/MouseWGBS_APCminTumorVsNormal/out_pileup/pileup.vcf.gz vs test/MouseWGBS_APCminTumorVsNormal/golden_pileup/pileup.vcf.gz
##compare test/MouseWGBS_APCminTumorVsNormal/out_pileup/pileup.vcf.stats vs test/MouseWGBS_APCminTumorVsNormal/golden_pileup/pileup.vcf.stats
-->

# single pileup

### pileup mouse WGBS Smadh3 chromosome 19 and chromosome M
<!---
$$$ mkdir -p test/Smadh3_chr19_chrM/out_pileup
@@ biscuit : biscuit-develop
@@ mouseSmadh3.bam : test/Smadh3_chr19_chrM/raw_bam/WGBS_Smadh3_chr19_chrM.bam
@@ mm10.fa : /home/wanding.zhou/references/mm10/mm10.fa
@@ mouseSmadh3.vcf : test/Smadh3_chr19_chrM/out_pileup/pileup.vcf
-->
```bash
$$ biscuit pileup -i mouseSmadh3.bam -r mm10.fa -g chr19:1-3100000 -o mouseSmadh3.vcf -q 3
```
<!---
##compare test/Smadh3_chr19_chrM/out_pileup/pileup.vcf vs test/Smadh3_chr19_chrM/golden_pileup/pileup.vcf
##compare test/Smadh3_chr19_chrM/out_pileup/pileup.vcf.stats vs test/Smadh3_chr19_chrM/golden_pileup/pileup.vcf.stats
-->

### pileup human TCGA WGBS NIC1254A46 chromosome 19 and chromosome M

<!---
$$$ mkdir -p test/NIC1254A46_chr19_chrM/out_pileup
@@ biscuit : biscuit-develop
@@ NIC1254A46.bam : test/NIC1254A46_chr19_chrM/raw_bam/NIC1254A46_chr19_chrM.bam
@@ hg19.fa : /home/wanding.zhou/references/hg19/hg19.fa
@@ NIC1254A46.vcf : test/NIC1254A46_chr19_chrM/out_pileup/pileup.vcf
-->
```bash
$$ biscuit pileup -i NIC1254A46.bam -r hg19.fa -o NIC1254A46.vcf -q 10
```
<!---
##compare test/NIC1254A46_chr19_chrM/out_pileup/pileup.vcf vs test/NIC1254A46_chr19_chrM/golden_pileup/pileup.vcf
##compare test/NIC1254A46_chr19_chrM/out_pileup/pileup.vcf.stats vs test/NIC1254A46_chr19_chrM/golden_pileup/pileup.vcf.stats
-->
