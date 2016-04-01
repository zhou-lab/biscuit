/**
 * Circulat Binary Segmentation
 * first perform ternary segmentation among blocks
 * then refine change points within blocks
 ***/

#include "wzcbs.h"

double *compute_partialsums(double *dat, int n, block_t *bk, psum_t *gpsum) {

  /**
   * ps[i] = ps[0] + ps[1] + ... + ps[i-1] is the partial sum before current
   * ps[0] === 0.0 is the change point before the first data
   * ps[n] is the change point after the last data
   **/
  double *ps = calloc(n+1, sizeof(double));

  gpsum->min_index = 0;
  gpsum->max_index = 0;
  gpsum->min = dat[0];
  gpsum->max = dat[0];

  int i; double p=0.0;
  for (i=0; i<n; ++i) {
    ps[i] = p;
    p += dat[i];
    if (p > gpsum->max) {
      gpsum->max = p;
      gpsum->max_index = i;
    }
    if (p < gpsum->min) {
      gpsum->min = p;
      gpsum->min_index = i;
    }

    block1_t *b = bk->a+i/bk->bsize;
    if (i == b->beg || p > b->psum.max) {
      b->psum.max = p;
      b->psum.max_index = i;
    }
    if (i == b->beg || p < b->psum.min) {
      b->psum.min = p;
      b->psum.min_index = i;
    }
  }
  ps[n] = p;

  return ps;
}

static inline double t_stats0(double rnjov1, double dpsum) {
  return rnjov1*dpsum*dpsum;
}

/* assert(alen >= 0) i.e., non-negative arc length */
static inline double rnjov1(int n, int alen) {
  return (double) n / (double) (alen*(n-alen));
}

static inline double t_stats(int n, int alen, double dpsum) {
  return t_stats0(rnjov1(n, alen), dpsum);
}

int compare_block_pairs(const void *_bp1, const void *_bp2) {
  if (((block_pair_t*) _bp1)->t > ((block_pair_t*) _bp2)->t) return 1;
  else if (((block_pair_t*) _bp1)->t < ((block_pair_t*) _bp2)->t) return -1;
  else return 0;
}

static void __debug_refine_tmax_fix_alen(int n, double *ps, block_pair_t *bp, int i2j) {
  int i;
  double t_max=-1; int i_max, j_max;
  for (i=MAX(bp->b->beg, bp->d->beg-i2j); i<=MIN(bp->b->end, bp->d->end-i2j); ++i) {
    int j = i+i2j;
    double t = t_stats(n, i2j+1, ps[i]-ps[j+1]);
    fprintf(stderr, "\t|| %d-%d: %1.3f", i, j ,t);
    if (t_max<0 || t_max < t) {
      t_max = t;
      i_max = i;
      j_max = j;
    }
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "i2j: %d || %d-%d : %1.3f\n", i2j, i_max, j_max, t_max);
}

/* [0....i-1][i....j][j+1....n-1] */
static void refine_tmax_fix_alen(int n, double *ps, block_pair_t *bp, int i2j, int chnpnts[2], double *t_gmax) {
  int i;
  for (i=MAX(bp->b->beg, bp->d->beg-i2j); i<=MIN(bp->b->end, bp->d->end-i2j); ++i) {
    int j = i+i2j;
    double t = t_stats(n, i2j+1, ps[i]-ps[j+1]);
    if (t > *t_gmax) {
      *t_gmax = t;
      chnpnts[0] = i;
      chnpnts[1] = j;
    }
  }
}

double *center_data(int *raw_dat, int n) {
  int64_t sum = 0;
  int i;
  for (i=0; i<n; ++i) sum += raw_dat[i];
  double mean = (double) sum / (double) n;
  double *dat = calloc(n, sizeof(double));
  for (i=0; i<n; ++i) dat[i] = (double) raw_dat[i] - mean;
  return dat;
}

/**
 * ternary segment (once)
 * number of change points can be 0,1,2
 * Assume data is centered! sum(dat) = 0
 * returned change point is identified by the position after
 * i.e., for change point [...c-1][c...] c is the recorded position */
void ternary_segmentation(int *raw_dat, int n, int *n_chnpnts, int chnpnts[2], double *t_gmax) {

  int verbose = 0;
  
  block_t *bk = init_blocks(n);
  psum_t gpsum;
  double *dat = center_data(raw_dat, n);
  double *ps = compute_partialsums(dat, n, bk, &gpsum);
  if (gpsum.min == gpsum.max) {
    *n_chnpnts = 0;
    goto ternary_segmentation_cleanup;
  }

  double t_gmax0 = t_stats(n, abs(gpsum.min_index-gpsum.max_index), gpsum.min-gpsum.max);

  /* identify block pairs that contains potential change points
   * this is the change point based on maximum partial sum difference */
  block_pair_v *block_pairs = init_block_pair_v(2);
  int bi, bj;
  for (bi=0; bi<bk->n; ++bi) {
    for (bj=bi; bj<bk->n; ++bj) {
      block1_t *b = bk->a+bi;
      block1_t *d = bk->a+bj;
      int alenhi = d->end - b->beg;
      int alenlo = bi==bj ? 1 : d->beg - b->end;
      double dpsum_bmd = b->psum.max - d->psum.min;
      double dpsum_dmb = d->psum.max - b->psum.min;
      double t_max = t_stats0((double) n/(double) MIN(alenhi*(n-alenhi), alenlo*(n-alenlo)), MAX(dpsum_bmd, dpsum_dmb));
      if (t_max >= t_gmax0) {
        block_pair_t *p = next_ref_block_pair_v(block_pairs);
        p->b = b;
        p->d = d;
        p->t_max = t_max;
        if (dpsum_bmd > dpsum_dmb) {
          p->alen = abs(b->psum.max_index - d->psum.min_index);
          p->t = t_stats(n, p->alen, dpsum_bmd);
        } else {
          p->alen = abs(d->psum.max_index - b->psum.min_index);
          p->t = t_stats(n, p->alen, dpsum_dmb);
        }
      }
    }
  }

  /* order by t-statistics */
  qsort(block_pairs->buffer, block_pairs->size, sizeof(block_pair_t), compare_block_pairs);

  *t_gmax = -1.0;
  /* refine change points inside the block pairs
   * this identifies the t-statistic max from partial-sum-difference max
   * to maximize t, we minimize alen*(n-alen), this is monotonous when alen is >n/2 and when alen is <n/2 */
  int bpi, i2j;
  for (bpi=0; bpi<block_pairs->size; ++bpi) {
    block_pair_t *bp = ref_block_pair_v(block_pairs, bpi);
    int alenhi = bp->d->end - bp->b->beg;
    int alenlo = bp->b==bp->d ? 0 : bp->d->beg-bp->b->end;

    if (bp->alen < n-1) {       /* avoid when two change points fuse */
      refine_tmax_fix_alen(n, ps, bp, bp->alen, chnpnts, t_gmax);
      if (verbose > 1) __debug_refine_tmax_fix_alen(n, ps, bp, bp->alen);
    }

    /* when alenlo < n/2, there is chance of minimizing alen*(n-alen) by decreasing alen */
    if (alenlo < n/2)
      for (i2j=alenlo; i2j<bp->alen; ++i2j) {
        refine_tmax_fix_alen(n, ps, bp, i2j, chnpnts, t_gmax);
        if (verbose > 1) __debug_refine_tmax_fix_alen(n, ps, bp, i2j);
      }

    /* when alenhi > n/2, there is chance of minimizing alen*(n-alen) by increasing alen */
    if (alenhi > n/2)
      for (i2j=alenhi; i2j>bp->alen; --i2j) {
        if (i2j >= n-1) continue; /* skip when two change points fuse */
        refine_tmax_fix_alen(n, ps, bp, i2j, chnpnts, t_gmax);
        if (verbose > 1) __debug_refine_tmax_fix_alen(n, ps, bp, i2j);
      }
  }

  /* determine the number of change points */
  if (chnpnts[0] == 0) {
    if (chnpnts[1] == n-1) {
      *n_chnpnts = 0;
    } else {
      /* [0...b][b+1...n-1]
       * chnpnts[0] is 0
       * chnpnts[1] is b and recorded as b+1 */
      *n_chnpnts = 1;
      chnpnts[0] = chnpnts[1]+1;
    }
  } else if (chnpnts[1] == n-1) {
    /* [0...a-1][a...n-1]
     * chnpnts[0] is a
     * chnpnts[1] is n-1 */
    *n_chnpnts = 1;
  } else {
    /* [0...a-1][a...b][b+1...n-1]
     * chnpnts[0] is a
     * chnpnts[1] is b and recorded as b+1 */
    *n_chnpnts = 2;
    chnpnts[1] += 1;
  }

 ternary_segmentation_cleanup:
  free(dat); free(ps); free(bk);
  return;
}


/* recursively segmentation, allocate segends
 * when no segmentation occurs, segends only have n
 * and n_segends == 1 */
int *recursive_segmentation(int *dat, int n, int *_n_segends) {

  int *segends = NULL;
  int n_segends = 0;

  int *ends = calloc(2, sizeof(int));
  ends[0] = 0; ends[1] = n;
  int k=1;                      /* k is always size(ends)-1 */
  int n_chnpnts; int chnpnts[2];
  double t_gmax;
  while (k>0) {
    ternary_segmentation(dat+ends[k-1], ends[k]-ends[k-1], &n_chnpnts, chnpnts, &t_gmax);
    if (n_chnpnts == 0) {
      /* prepend the segment end */
      segends = realloc(segends, sizeof(int)*(n_segends+1));
      if (n_segends>0) memmove(segends+1, segends, sizeof(int)*n_segends);
      segends[0] = ends[k];
      n_segends++;
      --k;
    } else if (n_chnpnts == 1) {
      ends = realloc(ends, sizeof(int)*(k+2));
      ends[k+1] = ends[k];
      ends[k] = chnpnts[0];
      k++;
    } else if (n_chnpnts == 2) {
      ends = realloc(ends, sizeof(int)*(k+3));
      ends[k+2] = ends[k];
      ends[k] = chnpnts[0];
      ends[k+1] = chnpnts[1];
      k += 2;
    }
  }
  free(ends);
  *_n_segends = n_segends;

  return segends;
}

int main(int argc, char *argv[])
{
  int n, i;
  int dat[100];
  int n_chnpnts;
  int chnpnts[2];
  double t_gmax;
  int kcn;
  int kcn1;

  /* n=22; */
  /* for (kcn=0; kcn<=n; ++kcn) { */
  /*   for (i=0; i<kcn; ++i) dat[i] = 4; */
  /*   for (i=kcn; i<n; ++i) dat[i] = 6; */
  /*   fprintf(stderr, "n: %d\n", n); */
  /*   ternary_segmentation(dat, n, &n_chnpnts, chnpnts, &t_gmax); */
  /*   fprintf(stderr, "nchnpnts: %d\n", n_chnpnts); */
  /*   fprintf(stderr, "chnpnts %d, %d\n", chnpnts[0], chnpnts[1]); */
  /*   fprintf(stderr, "expect: %d\n\n", kcn); */
  /*   /\* break; *\/ */
  /* } */

  /* n=22; */
  /* for (kcn=0; kcn<n-1; ++kcn) { */
  /*   for (kcn1=kcn+1; kcn1<n-1; ++kcn1) { */
  /*     /\* if (!(kcn==6 && kcn1==12)) continue; *\/ */
  /*     for (i=0; i<kcn; ++i) dat[i] = 4; */
  /*     for (i=kcn; i<kcn1; ++i) dat[i] = 6; */
  /*     for (i=kcn1; i<n; ++i) dat[i] = 4; */
  /*     fprintf(stderr, "n: %d\n",n); */
  /*     for (i=0; i<n; ++i) putchar('0'+dat[i]); */
  /*     putchar('\n'); */
  /*     ternary_segmentation(dat, n, &n_chnpnts, chnpnts, &t_gmax); */
  /*     fprintf(stderr, "nchnpnts: %d\n", n_chnpnts); */
  /*     fprintf(stderr, "chnpnts %d, %d\n", chnpnts[0], chnpnts[1]); */
  /*     fprintf(stderr, "expect: %d, %d\n\n", kcn, kcn1); */
  /*   } */
  /* } */

  /* recursive segmentation */
  /* n=22; */
  /* for (kcn=0; kcn<=n; ++kcn) { */
  /*   for (i=0; i<kcn; ++i) dat[i] = 4; */
  /*   for (i=kcn; i<n; ++i) dat[i] = 6; */
  /*   fprintf(stderr, "n: %d\n", n); */
  /*   int n_segends; */
  /*   int *segends = recursive_segmentation(dat, n, &n_segends); */
  /*   fprintf(stderr, "n_segends: %d\n", n_segends); */
  /*   fprintf(stderr, "segends: 0"); */
  /*   for (i=0; i<n_segends; ++i) { */
  /*     fprintf(stderr, ",%d", segends[i]); */
  /*   } */
  /*   putchar('\n'); */
  /*   fprintf(stderr, "expect: %d\n\n", kcn); */
  /*   free(segends); */
  /* } */

  n=22;
  for (kcn=0; kcn<n-1; ++kcn) {
    for (kcn1=kcn+1; kcn1<n-1; ++kcn1) {
      /* if (!(kcn==6 && kcn1==12)) continue; */
      for (i=0; i<kcn; ++i) dat[i] = 4;
      for (i=kcn; i<kcn1; ++i) dat[i] = 6;
      for (i=kcn1; i<n; ++i) dat[i] = 4;
      fprintf(stderr, "n: %d\n", n);
      int n_segends;
      int *segends = recursive_segmentation(dat, n, &n_segends);
      fprintf(stderr, "n_segends: %d\n", n_segends);
      fprintf(stderr, "segends: 0");
      for (i=0; i<n_segends; ++i) {
        fprintf(stderr, ",%d", segends[i]);
      }
      putchar('\n');
      fprintf(stderr, "expect: %d, %d\n\n", kcn, kcn1);
      free(segends);
    }
  }
  
  return 0;
}
