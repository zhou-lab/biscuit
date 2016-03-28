#ifndef WZCBS_H
#define WZCBS_H

/**
 * Circular Binary Segmentation 
 ***/

#include "wvec.h"
#include "math.h"

/* for grouping segments when there are too many segments */
typedef struct block_t {
  int beg;
  int end;
  int psmax_index;
  int psmin_index;
  int psmax;
  int psmin;
} block_t;

DEFINE_VECTOR(block_v, block_t)

typedef struct block_p {
  block_t *a, *b;
} block_p;

#endif /* WZCBS_H */
