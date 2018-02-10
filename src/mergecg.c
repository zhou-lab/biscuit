/* Merge C annd G beta values in CpG dinucleotide context
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2017 Wanding.Zhou@vai.org
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
 * 
**/

#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <zlib.h>
#include "wzbed.h"
#include "refcache.h"

typedef struct bed_data_meth_t {
  char ref;
  int nsamples;
  double *c_betas;
  double *g_betas;
  int *c_depts;
  int *g_depts;
} bed_data_meth_t;

void init_data_meth(bed1_t *b, void *aux_data) {
  (void) aux_data;
  b->data = calloc(1, sizeof(bed_data_meth_t));
}

void free_data_meth(void *bed_data) {
  free(((bed_data_meth_t*) bed_data)->c_betas);
  free(((bed_data_meth_t*) bed_data)->c_depts);
  free(((bed_data_meth_t*) bed_data)->g_betas);
  free(((bed_data_meth_t*) bed_data)->g_depts);
  free(bed_data);
}

void parse_data_meth(bed1_t *b, char **fields, int nfields) {

  int start;
  if (strcmp(fields[3],"C")==0 || strcmp(fields[3],"G")==0) start = 7;
  else start = 3;

  bed_data_meth_t *bd = (bed_data_meth_t*) b->data;
  if (bd->nsamples <= 0) { // first line
    bd->nsamples = (nfields-start)/2;
    if (bd->nsamples <= 0) wzfatal("No sample data identified.\n");
    bd->c_betas = calloc(bd->nsamples, sizeof(double));
    bd->c_depts = calloc(bd->nsamples, sizeof(int));
    bd->g_betas = calloc(bd->nsamples, sizeof(double));
    bd->g_depts = calloc(bd->nsamples, sizeof(int));
  } else if (bd->nsamples*2 + start != nfields) wzfatal("Malformed bed file.\n");
  
  // assume it's a C
  int i;
  for (i=0; i<bd->nsamples; ++i) {
    bd->c_betas[i] = atof(fields[start+2*i]);
    bd->c_depts[i] = atoi(fields[start+1+2*i]);
    memset(bd->g_betas, 0, bd->nsamples*sizeof(double));
    memset(bd->g_depts, 0, bd->nsamples*sizeof(int));
  }
}

static int usage() {

  fprintf(stderr, "\n");
  fprintf(stderr, "Merge C annd G beta values in CpG dinucleotide context\n");
  fprintf(stderr, "Usage: biscuit mergecg <ref.fa> <in.bed>\n");
  fprintf(stderr, "Input options:\n\n");
  fprintf(stderr, "     REF       fai-indexed fasta\n");
  fprintf(stderr, "     BED       sorted bed file, starting from column 4,5, we have alternating beta value and coverage for each sample. This is what output from biscuit vcf2bed without -e.\n");
  fprintf(stderr, "     -n        NOMe-seq mode, only merge C,G both in HCGD context\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");

  return 1;
}

void format_output(bed1_t *p, char *chrm, char base_before, char base_after) {
  // skip reporting if all sample has zero coverage
  int i; int allzero_cov = 1;
  bed_data_meth_t *pd = (bed_data_meth_t*) p->data;
  for (i=0; i<pd->nsamples; ++i) {
    if (pd->c_depts[i] + pd->g_depts[i] > 0) {
      allzero_cov = 0; break;
    }
  }
  if (allzero_cov) return;

  // if p is in CpG contxt, adjust the coordinates
  if (pd->ref == 'C' && base_after == 'G') p->end++;
  else if (pd->ref == 'G' && base_before == 'C') p->beg--;

  printf("%s\t%"PRId64"\t%"PRId64, chrm, p->beg, p->end);
  for (i=0; i<pd->nsamples; ++i) {
    int cov = pd->c_depts[i] + pd->g_depts[i];
    if (cov == 0) {
      fputs("\t.\t0", stdout);
    } else {
      printf("\t%1.3f\t%d", (pd->c_betas[i]*pd->c_depts[i] + pd->g_betas[i]*pd->g_depts[i]) / (double) cov, cov);
    }
    if (pd->c_depts[i] == 0) {
      fputs("\tC:.:0", stdout);
    } else {
      printf("\tC:%1.3f:%d", pd->c_betas[i], pd->c_depts[i]);
    }
    if (pd->g_depts[i] == 0) {
      fputs(",G:.:0", stdout);
    } else {
      printf(",G:%1.3f:%d", pd->g_betas[i], pd->g_depts[i]);
    }
  }
  putchar('\n');
}

int main_mergecg(int argc, char *argv[]) {

  int c;
  int nome_mode = 0;
  while ((c = getopt(argc, argv, "nh"))>=0) {
    switch (c) {
      case 'n': nome_mode=1; break;
      case 'h': return usage(); break;
      default: wzfatal("Unrecognized option: %c.\n", c); break;
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply reference file and sorted bed files.\n"); 
  }

  refcache_t *rc = init_refcache(argv[optind++], 10000, 10000);
  bed_file_t *bed = init_bed_file(argv[optind++]);
  bed1_t *b = init_bed1(init_data_meth, NULL);
  bed1_t *p = init_bed1(init_data_meth, NULL);
  b->tid = -1; p->tid = -1;
  bed_data_meth_t *bd, *pd;
  char p_base_before='N', p_base_after='N';
  char b_base_before='N', b_base_after='N';
  while (bed_read1(bed, b, parse_data_meth)) {
    refcache_fetch_chrm(rc, tid2name(bed->targets, b->tid));
    bd = (bed_data_meth_t*) b->data;
    bd->ref = refcache_getbase_upcase(rc, b->end);
    b_base_before = refcache_getbase_upcase(rc, b->end-1);
    b_base_after = refcache_getbase_upcase(rc, b->end+1);
    if (bd->ref == 'G') { // correct based on the actual reference base
      memcpy(bd->g_betas, bd->c_betas, bd->nsamples*sizeof(double));
      memcpy(bd->g_depts, bd->c_depts, bd->nsamples*sizeof(int));
      memset(bd->c_betas, 0, bd->nsamples*sizeof(double));
      memset(bd->c_depts, 0, bd->nsamples*sizeof(int));
    }

    // merge C, G
    pd = (bed_data_meth_t*) p->data;
    if (b->tid == p->tid && 
        b->beg == p->beg+1 && b->end == p->end+1 && 
        bd->ref == 'G' && pd->ref == 'C' &&
        (!nome_mode ||
         (p_base_before != 'G' && b_base_after != 'C'))) {
      if (pd->nsamples != bd->nsamples) wzfatal("Missing sample at %s:%u-%u.\n", rc->chrm, b->beg, b->end);
      memcpy(pd->g_betas, bd->g_betas, bd->nsamples*sizeof(double));
      memcpy(pd->g_depts, bd->g_depts, bd->nsamples*sizeof(int));
      b->tid = -1; // merged, no more reporting
    }

    // report p
    if (p->tid >= 0) format_output(p, tid2name(bed->targets, p->tid), p_base_before, p_base_after);
    bed1_t *tmp = p; p = b; b = tmp;
    p_base_before = b_base_before; p_base_after = b_base_after;
  }

  // report last p
  if (p->tid >= 0) format_output(p, tid2name(bed->targets, p->tid), p_base_before, p_base_after);
  free_bed1(b, free_data_meth);
  free_bed1(p, free_data_meth);
  free_refcache(rc);
  free_bed_file(bed);
  
  return 0;
}
