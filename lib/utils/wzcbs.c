/**
 * Circulat Binary Segmentation
 ***/

int* compute_partialsums(int *dat, int n, block_t *bk, psum_t *gpsum) {

  int ps = calloc(n, sizeof(int));
  
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

double _t_stats(double rnjov1, int dpsum) {
  return rnjov1*(dpsum**2);
}

/* assert(alen >= 0) i.e., non-negative arc length */
double rnjov1(int n, int alen) {
  return (double) n / (double) alen / (double) (n-alen);
}

double t_stats(int n, int alen, int psi, int psj) {
  return _t_stats(rnjov1(n, alen), psi-psj);
}

int compare_block_pairs(void *_bp1, void *_bp2) {
  if (((block_pair_t*) _bp1)->b->t_stats > ((block_pair_t*) _bp2)->d->t_stats) return 1;
  else if (((block_pair_t*) _bp1)->b->t_stats < ((block_pair_t*) _bp2)->d->t_stats) return -1;
  else return 0;
}

void refine_tmax_fix_alen(int *ps, block_pair_t *bp, int i2j, int chnpnts[2], double *t_stats_gmax) {
  int i;
  for (i=max(bp->b->beg, bp->d->beg-i2j); i<=min(bp->b->end, bp->d->end-i2j); ++i) {
    int j = i+i2j;
    double t_stats = t_stats(n, i2j, ps[i], ps[j]);
    if (t_stats > *t_stats_gmax) {
      *t_stats_gmax = t_stats;
      chnpnts[0] = i;
      chnpnts[1] = j;
    }
  }
}

/* ternary segment (once)
 * number of change points can be 0,1,2 */
void seg_once(int *dat, int n, int *n_chnpnts, int chnpnts[2], double *t_stats_gmax) {

  block_t *bk = init_blocks(n);
  psum_t gpsum;
  int *ps = compute_partialsums(dat, n, bk, &gpsum);
  double t_stats_gmax = t_stats(n, abs(gpsum.min_index-gpsum.max_index), gpsum.min, gpsum.max);

  /* identify block pairs that contains potential change points
   * this is the change point based on maximum partial sum difference */
  block_pair_v *block_pairs = init_block_pair_v(2);
  int bi, bj;
  for (bi=0; bi<nb; ++bi) {
    for (bj=bi; bj<nb; ++bj) {
      block1_t *b = bk->a+bi;
      block1_t *d = bk->a+bj;
      int alenhi = d->end - b->beg;
      int alenlo = d->beg - b->end;
      int dpsum_bmd = b->psum.max - d->psum.min;
      int dpsum_dmb = d->psum.max - b->psum.min;
      double t_stats_max = _t_stats((double) n/(double) MIN(alenhi*(n-alenhi), alenlo*(n-alenlo)),
                                    MAX(dpsum_bmd, dpsum_dmb));
      if (t_stats_max >= t_stats_gmax) {
        block_pair_t *p = next_ref_block_pair_v(block_pairs);
        p->b = b;
        p->d = d;
        p->t_stats_max = t_stats_max;
        if (dpsum_bmd > dpsum_dmb) {
          p->alen = abs(b->psum.max_index - d->psum.min_index);
          p->t_stats = t_stats(n, p->alen, dpsum_bmd);
        } else {
          p->alen = abs(d->psum.max_index - b->psum.min_index);
          p->t_stats = t_stats(n, p->alen, dpsum_dmb);
        }
      }
    }
  }

  /* order by t-statistics */
  qsort(block_pairs->buffer, block_pairs->size, sizeof(block_pair_t), compare_block_pairs);

  /* refine change points inside the block pairs
   * this identifies the t-statistic max from partial-sum-difference max
   * to maximize t, we minimize alen*(n-alen), this is monotonous when alen is >n/2 and when alen is <n/2 */
  for (i=0; i<block_pairs->size; ++i) {
    block_pair_t *bp = ref_block_pair_v(block_pairs, i);
    alenhi = bp->d->end - bp->b->beg;
    alenlo = bp->d->beg - bp->b->end;

    /* when alenlo < n/2, there is chance of minimizing alen*(n-alen) by decreasing alen */
    if (alenlo < n/2)
      for (i2j=alenlo; i2j<alenmax; ++i2j)
        refine_tmax_fix_alen(ps, bp, i2j, chnpnts, t_stats_gmax);

    /* when alenhi > n/2, there is chance of minimizing alen*(n-alen) by increasing alen */
    if (alenhi > n/2)
      for (i2j=alenhi; i2j>alenmax; --i2j)
        refine_tmax_fix_alen(ps, bp, i2j, chnpnts, t_stats_gmax);
  }

  /* determine the number of change points */
  if (chnpnts[0] == 0) {
    if (chnpnts[1] == n-1) {
      *nchnpnts = 0;
    } else {
      *nchnpnts = 1;
      chnpnts[0] = chnpnts[1];
    }
  } else if (chnpnts[1] == n-1) {
    *nchnpnts = 1;
  } else {
    *nchnpnts = 2;
  }
}


/* recursively segmentation, allocate segends
 * when no segmentation occurs, segends only have n
 * and n_segends == 1 */
int seg_recursive(int *dat, int n, int *segends) {

  segends = NULL;
  int n_segends = 0;

  int *ends = calloc(2, sizeof(int));
  ends[0] = 0; ends[1] = n;
  int k=1;                      /* k is always size(ends)-1 */
  int n_chnpnts; int chnpnts[2];
  while (k>0) {
    seg_once(dat[ends[k-1]], ends[k-1]-ends[k], &n_chnpnts, chnpnts);
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
