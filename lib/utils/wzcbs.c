/**
 * Circulat Binary Segmentation
 ***/

void compute_partialsums() {

}

/* ternary segment (once)
 * number of change points can be 0,1,2 */
void seg_once(int *dat, int n, int *n_chnpnts, int chnpnts[2]) {

  /* set blocks */
  int nb = n>50?int(sqrt((double)n)):1;
  int bsize = (n-1)/nb+1;
  block_t *blocks = calloc(nb, sizeof(block_t));
  int bi;
  block_t *b;
  for (bi=0; bi<nb; ++bi) {
    b = blocks+bi;
    b->beg = bi*bsize;
    b->end = (bi+1)*bsize-1;
  }
  blocks[nb-1].end = n-1;

  int ps = calloc(n, sizeof(int)); /* partial sum */
  int psmin_index = 0;
  int psmax_index = 0;
  int psmin = dat[0];
  int psmax = dat[0];

  int i; int p=0;
  for (i=0; i<n; ++i) {
    p += dat[i];
    ps[i] = p;
    if (p > psmax) {
      psmax = p;
      psmax_index = i;
    }
    if (p < psmin) {
      psmin = p;
      psmin_index = i;
    }

    bi = i/bsize;
    b = blocks+bi;
    if (i == b->beg || p > b->psmax) {
      b->psmax = p;
      b->psmax_index = i;
    }
    if (i == b->beg || p < b->psmin) {
      b->psmin = p;
      b->psmin_index = i;
    }
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
