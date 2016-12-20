/**
 * main.c
 * The MIT License (MIT)
 * Copyright (c) 2016 Wanding.Zhou@vai.org
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
/* #include "utils.h" */
#include "biscuit.h"

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
int main_markdup(int argc, char *argv[]);
int main_ndr(int argc, char *argv[]);
int main_vcf2bed(int argc, char *argv[]);
int main_epiread(int argc, char *argv[]);
int main_asm(int argc, char *argv[]);

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: biscuit (bisulfite data tool kit)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Wanding Zhou <wanding.zhou@vai.org>\n\n");
  fprintf(stderr, "Usage:   biscuit <command> [options]\n\n");
  fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
  fprintf(stderr, "         align         align bisulfite treated short reads using adapted BWA-mem algorithm\n");
  fprintf(stderr, "         markdup       mark duplicates on the same bisulfite strand\n");
  fprintf(stderr, "         pileup        pileup cytosine and mutations, estimate bisulfite conversion rate and meth-level averages.\n");
  fprintf(stderr, "         ndr           call nucleosome depletion region (NDR) from NOMe-seq\n");
  fprintf(stderr, "         vcf2bed       convert VCF to bed graph\n");
  fprintf(stderr, "         epiread       convert bam to epiread format\n");
  fprintf(stderr, "         asm           test allele specific methylation\n");
  /* fprint */
  /* fprintf(stderr, "         hemi          find hemi-methylated region.\n"); */
  fprintf(stderr, "\n");

  return 1;
}

int main(int argc, char *argv[])
{
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
  else if (strcmp(argv[1], "markdup") == 0) ret = main_markdup(argc-1, argv+1);
  else if (strcmp(argv[1], "ndr") == 0) ret = main_ndr(argc-1, argv+1);
  else if (strcmp(argv[1], "vcf2bed") == 0) ret = main_vcf2bed(argc-1, argv+1);
  else if (strcmp(argv[1], "epiread") == 0) ret = main_epiread(argc-1, argv+1);
  else if (strcmp(argv[1], "asm") == 0) ret = main_asm(argc-1, argv+1);
  else {
    fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
    return 1;
  }

  fflush(stdout);               /* not enough for remote file systems */
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
