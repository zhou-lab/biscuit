#ifndef WZCBS_H
#define WZCBS_H

/**
 * Circular Binary Segmentation 
 ***/

#include "wvec.h"
#include "math.h"

#define MAX(a,b)                \
  ({ __typeof__ (a) _a = (a);   \
    __typeof__ (b) _b = (b);    \
    _a > _b ? _a : _b; })

#define MIN(a,b)                \
  ({ __typeof__ (a) _a = (a);   \
    __typeof__ (b) _b = (b);    \
    _a > _b ? _b : _a; })

typedef struct psum_t {
  int max;
  int min;
  int max_index;
  int min_index;
} psum_t;

/* for grouping segments when there are too many segments */
typedef struct block1_t {
  int beg;
  int end;
  psum_t psum;
} block1_t;


typedef struct block_pair_t {
  block1_t *b, *d;
  double t, t_max;
  int alen;                  /* arc length */
} block_pair_t;

DEFINE_VECTOR(block_pair_v, block_pair_t)

typedef struct block_t {
  block1_t *a;
  int n;
  int bsize;
} block_t;

static inline block_t *init_blocks(int n) {

  block_t *bk = calloc(1,sizeof(block_t));

  bk->n = n>50 ? (int)sqrt((double)n) : 1;
  bk->bsize = (n-1) / bk->n+1;
  bk->a = calloc(bk->n, sizeof(block1_t));
  int bi;
  for (bi=0; bi<bk->n; ++bi) {
    block1_t *b = bk->a+bi;
    b->beg = bi*bk->bsize;
    b->end = (bi+1)*bk->bsize-1;
  }
  bk->a[bk->n-1].end = n-1;

  return bk;
}

#endif /* WZCBS_H */
