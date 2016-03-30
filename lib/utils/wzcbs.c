/**
 * Circulat Binary Segmentation
 * first perform ternary segmentation among blocks
 * then refine change points within blocks
 ***/

#include "wzcbs.h"

int *compute_partialsums(int *dat, int n, block_t *bk, psum_t *gpsum) {

  int *ps = calloc(n, sizeof(int));
  
  gpsum->min_index = 0;
  gpsum->max_index = 0;
  gpsum->min = dat[0];
  gpsum->max = dat[0];
  
  int i; int p=0;
  for (i=0; i<n; ++i) {
    p += dat[i];
    ps[i] = p;
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

  return ps;
}

static inline double t_stats0(double rnjov1, int dpsum) {
  return rnjov1*dpsum*dpsum;
}

/* assert(alen >= 0) i.e., non-negative arc length */
static inline double rnjov1(int n, int alen) {
  return (double) n / (double) alen / (double) (n-alen);
}

static inline double t_stats(int n, int alen, int dpsum) {
  return t_stats0(rnjov1(n, alen), dpsum);
}

int compare_block_pairs(const void *_bp1, const void *_bp2) {
  if (((block_pair_t*) _bp1)->t > ((block_pair_t*) _bp2)->t) return 1;
  else if (((block_pair_t*) _bp1)->t < ((block_pair_t*) _bp2)->t) return -1;
  else return 0;
}

void refine_tmax_fix_alen(int n, int *ps, block_pair_t *bp, int i2j, int chnpnts[2], double *t_gmax) {
  int i;
  for (i=MAX(bp->b->beg, bp->d->beg-i2j); i<=MIN(bp->b->end, bp->d->end-i2j); ++i) {
    int j = i+i2j;
    double t = t_stats(n, i2j, ps[i]-ps[j]);
    if (t > *t_gmax) {
      *t_gmax = t;
      chnpnts[0] = i;
      chnpnts[1] = j;
    }
  }
}

/* ternary segment (once)
 * number of change points can be 0,1,2
 * Assume data is centered! sum(dat) = 0 */
void ternary_segmentation(int *dat, int n, int *n_chnpnts, int chnpnts[2], double *t_gmax) {

  block_t *bk = init_blocks(n);
  psum_t gpsum;
  int *ps = compute_partialsums(dat, n, bk, &gpsum);
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
      int dpsum_bmd = b->psum.max - d->psum.min;
      int dpsum_dmb = d->psum.max - b->psum.min;
      double t_max = t_stats0((double) n/(double) MIN(alenhi*(n-alenhi), alenlo*(n-alenlo)),
                                    MAX(dpsum_bmd, dpsum_dmb));
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

  /* refine change points inside the block pairs
   * this identifies the t-statistic max from partial-sum-difference max
   * to maximize t, we minimize alen*(n-alen), this is monotonous when alen is >n/2 and when alen is <n/2 */
  int bpi, i2j;
  for (bpi=0; bpi<block_pairs->size; ++bpi) {
    block_pair_t *bp = ref_block_pair_v(block_pairs, bpi);
    int alenhi = bp->d->end - bp->b->beg;
    int alenlo = bp->b==bp->d ? 1 : bp->d->beg-bp->b->end;

    /* when alenlo < n/2, there is chance of minimizing alen*(n-alen) by decreasing alen */
    if (alenlo < n/2)
      for (i2j=alenlo; i2j<bp->alen; ++i2j)
        refine_tmax_fix_alen(n, ps, bp, i2j, chnpnts, t_gmax);

    /* when alenhi > n/2, there is chance of minimizing alen*(n-alen) by increasing alen */
    if (alenhi > n/2)
      for (i2j=alenhi; i2j>bp->alen; --i2j)
        refine_tmax_fix_alen(n, ps, bp, i2j, chnpnts, t_gmax);
  }

  /* determine the number of change points */
  if (chnpnts[0] == 0) {
    if (chnpnts[1] == n-1) {
      *n_chnpnts = 0;
    } else {
      *n_chnpnts = 1;
      chnpnts[0] = chnpnts[1];
    }
  } else if (chnpnts[1] == n-1) {
    *n_chnpnts = 1;
  } else {
    *n_chnpnts = 2;
  }
}


/* recursively segmentation, allocate segends
 * when no segmentation occurs, segends only have n
 * and n_segends == 1 */
int recursive_segmentation(int *dat, int n, int *segends) {

  segends = NULL;
  int n_segends = 0;

  int *ends = calloc(2, sizeof(int));
  ends[0] = 0; ends[1] = n;
  int k=1;                      /* k is always size(ends)-1 */
  int n_chnpnts; int chnpnts[2];
  double t_gmax;
  while (k>0) {
    ternary_segmentation(dat+ends[k-1], ends[k-1]-ends[k], &n_chnpnts, chnpnts, &t_gmax);
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

  return n_segends;
}

int main(int argc, char *argv[])
{
  /* int dat[20] = {4,4,4,4,4,4,4,4,4,4, */
  /*                6,6,6,6,6,6,6,6,6,6}; */

  int dat[20] = {1,1,1,1,1,1,1,1,1,1,
                 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  int n_chnpnts;
  int chnpnts[2];
  double t_gmax;
  ternary_segmentation(dat, 20, &n_chnpnts, chnpnts, &t_gmax);
  
  /* recursive_segmentation() */
  return 0;
}
