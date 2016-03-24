/**
 * Circulat Binary Segmentation
 ***/

/* ternary segment (once)
 * number of change points can be 0,1,2 */
void seg_once(int *dat, int n, int *n_chnpnts, int chnpnts[2]) {
  
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
