#include "stats.h"
#include "wztsv.h"
#include "encode.h"

/* typedef struct { */
/*   int verbose; */
/*   /\* for emission probability *\/ */
/*   int a; */
/*   int b; */
/*   /\* for transition probability *\/ */
/*   double p; */
/*   double q; */
/* } conf_t; */

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
    fprintf(stdout, "%s\t%d\t%d\t%c/%c\t%c/%c\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n",
            chrm, snp_loc, cg_loc,
            nt256int8_to_nt256char_table[smax[0]],
            nt256int8_to_nt256char_table[smax[1]],
            nt256int8_to_nt256char_table[cmax[0]],
            nt256int8_to_nt256char_table[cmax[1]],
            cross[smax[0]*5+cmax[0]],
            cross[smax[0]*5+cmax[1]],
            cross[smax[1]*5+cmax[0]],
            cross[smax[1]*5+cmax[1]],
            left, right, two);
  }
}

int main_asm(int argc, char *argv[]) {

  /* conf_t conf = {}; */

  if (optind >= argc) {
    fprintf(stderr, "Please provide input epiread.\n");
    fflush(stderr);
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
  
  /* int c, i; int collapse=0; char *out_fn=0; */
  /* while ((c = getopt(argc, argv, "i:")) >= 0) { */
  /*   switch (c) { */
  /*   case 'V': conf.verbose = atoi(optarg); break; */
  /*   case 'b': in->bed = optarg; break; */
  /*   case 'o': out_fn = optarg; break; */
  /*   case 'A': conf.a = atoi(optarg); break; */
  /*   case 'B': conf.b = atoi(optarg); break; */
  /*   case 'P': conf.p = atof(optarg); break; */
  /*   case 'Q': conf.q = atof(optarg); break; */
  /*   case 'c': collapse = 1; break; */
  /*   case 'h': { */
  /*     fprintf(stderr, "\n"); */
  /*     fprintf(stderr, "Usage: biscuit asm [options] pairwise_epiread \n"); */
  /*     fprintf(stderr, "Input options:\n"); */
  /*     fprintf(stderr, "     -b FILE   bed file of GpC retention, coordinate-sorted\n"); */
  /*     fprintf(stderr, "     -c        collapse into regions\n"); */
  /*     fprintf(stderr, "     -A INT    parameter a in beta binomial emission [%d]\n", conf.a); */
  /*     fprintf(stderr, "     -B INT    parameter b in beta binomial emission [%d]\n", conf.b); */
  /*     fprintf(stderr, "     -P FLOAT  parameter p in transition probability [%1.2f]\n", conf.p); */
  /*     fprintf(stderr, "     -Q FLOAT  parameter q in transition probability [%1.2f]\n", conf.q); */
  /*     fprintf(stderr, "     -o FILE   output file\n"); */
  /*     fprintf(stderr, "     -V INT    verbose level [%d].\n", conf.verbose); */
  /*     fprintf(stderr, "     -h        this help.\n"); */
  /*     fprintf(stderr, "\n"); */
  /*     return 1; */
  /*   } */
  /*   default: */
  /*     fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c); */
  /*     fflush(stderr); */
  /*     exit(1); */
  /*     break; */
  /*   } */
  /* } */

  return 0;
}
