#ifndef WZCBS_H
#define WZCBS_H

/**
 * Circular Binary Segmentation 
 ***/

#include "wvec.h"

/* for grouping segments when there are too many segments */
typedef struct block_t {
  ;
} block_t;


typedef struct block_p {
  block_t *a, *b;
} block_p;

#endif /* WZCBS_H */
