/**
 * main.c
 * The MIT License (MIT)
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
**/

#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "kstring.h"
#include "biscuit.h"
#include "wzmisc.h"
//#include "../lib/htslib/version.h"

static double cputime() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + tp.tv_usec * 1e-6;
}

int main_biscuit_index(int argc, char *argv[]);
int main_align(int argc, char *argv[]);
int main_pileup(int argc, char *argv[]);
/* int main_ndr(int argc, char *argv[]); */
int main_vcf2bed(int argc, char *argv[]);
int main_epiread(int argc, char *argv[]);
int main_asm(int argc, char *argv[]);
int main_tview(int argc, char *argv[]);
int main_bsstrand(int argc, char *argv[]);
int main_cinread(int argc, char *argv[]);
int main_bsconv(int argc, char *argv[]);
int main_mergecg(int argc, char *argv[]);
int main_rectangle(int argc, char *argv[]);
int main_qc(int argc, char *argv[]);

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: BISCUIT (BISulfite-seq CUI Toolkit)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Jacob Morrison <jacob.morrison@vai.org>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: biscuit <command> [options]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Command:\n");
  fprintf(stderr, " -- Read mapping\n");
  fprintf(stderr, "    index        Index reference genome sequences in the FASTA format\n");
  fprintf(stderr, "    align        Align bisulfite treated short reads using adapted BWA-mem\n");
  fprintf(stderr, "                     algorithm\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " -- BAM operation\n");
  fprintf(stderr, "    tview        Text alignment viewer with bisulfite coloring\n");
  fprintf(stderr, "    bsstrand     Validate/correct bisulfite conversion strand label (YD tag)\n");
  fprintf(stderr, "    bsconv       Summarize/filter reads by bisulfite conversion (ZN tag)\n");
  fprintf(stderr, "    cinread      Print cytosine-read pair in a long form\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " -- Base summary\n");
  fprintf(stderr, "    pileup       Pileup cytosine and mutations\n");
  fprintf(stderr, "    vcf2bed      Convert VCF to BED file\n");
  fprintf(stderr, "    mergecg      Merge C and G in CpG context\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " -- Epireads\n");
  fprintf(stderr, "    epiread      Convert BAM to epiread format\n");
  fprintf(stderr, "    rectangle    Convert epiread to rectangle format\n");
  fprintf(stderr, "    asm          Test allele specific methylation\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " -- Other\n");
  fprintf(stderr, "    version      Print BISCUIT and library versions\n");
  fprintf(stderr, "    qc           Generate QC files to hook to MultiQC\n");
  /* fprintf(stderr, "    ndr          Call nucleosome depletion region (NDR) from NOMe-seq\n"); */
  fprintf(stderr, "\n");

  return 1;
}

int main(int argc, char *argv[]) {
  extern char *bwa_pg;
  int i, ret;
  double t_real;
  kstring_t pg = {0,0,0};
  t_real = realtime();
  ksprintf(&pg, "@PG\tID:biscuit\tPN:biscuit\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
  for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
  bwa_pg = pg.s;
  if (argc < 2) return usage();
  if (strcmp(argv[1], "index") == 0) ret = main_biscuit_index(argc-1, argv+1);
  else if (strcmp(argv[1], "align") == 0) ret = main_align(argc-1, argv+1);
  else if (strcmp(argv[1], "pileup") == 0) ret = main_pileup(argc-1, argv+1);
  /* else if (strcmp(argv[1], "ndr") == 0) ret = main_ndr(argc-1, argv+1); */
  else if (strcmp(argv[1], "vcf2bed") == 0) ret = main_vcf2bed(argc-1, argv+1);
  else if (strcmp(argv[1], "epiread") == 0) ret = main_epiread(argc-1, argv+1);
  else if (strcmp(argv[1], "asm") == 0) ret = main_asm(argc-1, argv+1);
  else if (strcmp(argv[1], "tview") == 0) ret = main_tview(argc-1, argv+1);
  else if (strcmp(argv[1], "bsstrand") == 0) ret = main_bsstrand(argc-1, argv+1);
  else if (strcmp(argv[1], "cinread") == 0) ret = main_cinread(argc-1, argv+1);
  else if (strcmp(argv[1], "mergecg") == 0) ret = main_mergecg(argc-1, argv+1);
  else if (strcmp(argv[1], "bsconv") == 0) ret = main_bsconv(argc-1, argv+1);
  else if (strcmp(argv[1], "rectangle") == 0) ret = main_rectangle(argc-1, argv+1);
  else if (strcmp(argv[1], "qc") == 0) ret = main_qc(argc-1, argv+1);
  else if (strcmp(argv[1], "version") == 0) {
      fprintf(stderr, "BISCUIT Version: %s\n\n", PACKAGE_VERSION);
      fprintf(stderr, "Using:\n");
      fprintf(stderr, "\thtslib version: zwdzwd/htslib at commit f29fa32\n");
      fprintf(stderr, "\tklib   version: zwdzwd/klib   at commit ca862f8\n");
      fprintf(stderr, "\tsgsl   version: zwdzwd/sgsl   at commit a0ddc77\n");
      fprintf(stderr, "\tutils  version: zwdzwd/utils  at commit 332459b\n");
      fprintf(stderr, "\nLibraries found at: https://github.com/huishenlab/biscuit/tree/master/lib\n");
      return 0;
  }
  else {
    usage();
    wzfatal("Unrecognized subcommand: %s\n", argv[1]);
  }

  fflush(stdout);             /* not enough for remote file systems */
  fclose(stdout);
  if (ret == 0) {
    fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (i = 0; i < argc; ++i)
      fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime()-t_real, cputime());
  }
  free(bwa_pg);

  return ret;
}
