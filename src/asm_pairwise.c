/* extract allele-specific methylation, pairwise
 * 
 * The MIT License (MIT)
 *
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
 * 
**/

#include "stats.h"
#include "wztsv.h"
#include "encode.h"
#include "gsl/gsl_cdf.h"

int *row_sums(int *matrix, int nr, int nc) {
  int *sums = calloc(nr, sizeof(int));
  int i, j;
  for (i=0; i<nr; ++i)
    for (j=0; j<nc; ++j)
      sums[i] += matrix[i*nc + j];
  return sums;
}

int *col_sums(int *matrix, int nr, int nc) {
  int *sums = calloc(nc, sizeof(int));
  int i, j;
  for (i=0; i<nr; ++i)
    for (j=0; j<nc; ++j)
      sums[j] += matrix[i*nc + j];
  return sums;
}

void max2(int *v, int n, int maxes[2]) {
  maxes[0] = 0; maxes[1] = 1;
  int i;
  for (i=2; i<n; ++i) {
    int minmax = v[maxes[0]] < v[maxes[1]] ? 0 : 1;
    if (v[i] >= v[maxes[minmax]])
      maxes[minmax] = i;
  }
}

/** fisher's exact test and chi-square test */
void test_asm(int *cross, char *chrm, int snp_loc, int cg_loc) {
  int *rs = row_sums(cross, 5, 5);
  int smax[2];
  max2(rs, 5, smax);
  int *cs = col_sums(cross, 5, 5);
  int cmax[2];
  max2(cs, 5, cmax);
  if (rs[smax[0]] > 0 && rs[smax[1]] > 0 &&
      cs[cmax[0]] > 0 && cs[cmax[1]] > 0) {
    double left, right, two;
    fisher_exact(cross[smax[0]*5+cmax[0]],
                 cross[smax[0]*5+cmax[1]],
                 cross[smax[1]*5+cmax[0]],
                 cross[smax[1]*5+cmax[1]],
                 &left, &right, &two);

    double pchisq = gsl_cdf_chisq_Q(two_by_two_chisq(cross[smax[0]*5+cmax[0]],
                                                     cross[smax[0]*5+cmax[1]],
                                                     cross[smax[1]*5+cmax[0]],
                                                     cross[smax[1]*5+cmax[1]]), 2);
    if (snp_loc != cg_loc &&
        nt256int8_to_nt256char_table[cmax[0]] != 'N' &&
        nt256int8_to_nt256char_table[cmax[1]] != 'N') {
      fprintf(stdout, "%s\t%d\t%d\t%c/%c\t%c/%c\t%d\t%d\t%d\t%d\t%f\t%f\n",
              chrm, snp_loc, cg_loc,
              nt256int8_to_nt256char_table[smax[0]],
              nt256int8_to_nt256char_table[smax[1]],
              nt256int8_to_nt256char_table[cmax[0]],
              nt256int8_to_nt256char_table[cmax[1]],
              cross[smax[0]*5+cmax[0]],
              cross[smax[0]*5+cmax[1]],
              cross[smax[1]*5+cmax[0]],
              cross[smax[1]*5+cmax[1]], two, pchisq);
    }
  }
}

static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit asm [options] <in.epiread>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -h    This help\n");
    fprintf(stderr, "\n");
}

int main_asm(int argc, char *argv[]) {

    int c;
    if (argc<2) { usage(); return 1; }
    while ((c = getopt(argc, argv, ":h")) >= 0) {
        switch(c) {
            case 'h': usage(); return 1;
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt);
            default: usage(); return 1;
        }
    }

  if (optind >= argc) {
    usage();
    wzfatal("Missing in.epiread");
  }

  char *input_fn = argv[optind];
  tsv_t *in = tsv_open(input_fn);

  char *chrm = NULL;
  int snp_loc=-1, cg_loc=-1, _snp_loc, _cg_loc;
  int cross[25];          /* row: snp, column: cg|hcg|gch */
  memset(cross, 0, 25*sizeof(int));
  while(tsv_read(in)) {
    if (in->n < 5) continue;    /* incomplete line */

    _snp_loc = atoi(in->fields[1]);
    _cg_loc = atoi(in->fields[2]);
    if (chrm == NULL || cg_loc != _cg_loc || snp_loc != _snp_loc ||
        strcmp(chrm, in->fields[0]) != 0) {
      if (chrm != NULL)
        test_asm(cross, chrm, snp_loc, cg_loc);
      chrm = realloc(chrm, strlen(in->fields[0])+1);
      strcpy(chrm, in->fields[0]);
      cg_loc = _cg_loc;
      snp_loc = _snp_loc;
      memset(cross, 0, 25*sizeof(int));
    }
    
    assert(strlen(in->fields[3]) == 1);
    assert(strlen(in->fields[4]) == 1);
    unsigned char _snp = in->fields[3][0];
    unsigned char _cg = in->fields[4][0];
    cross[nt256char_to_nt256int8_table[_snp]*5 +
          nt256char_to_nt256int8_table[_cg]]++;
  }
  if (chrm != NULL)
    test_asm(cross, chrm, snp_loc, cg_loc);

  free(chrm);
  tsv_close(in);
 
  return 0;
}



