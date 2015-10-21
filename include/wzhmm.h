#ifndef _WZHMM_H
#define _WZHMM_H

/**************************************
 * discrete state Markov chain (DSMC) *
 **************************************/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <zlib.h>
#include <inttypes.h>

typedef struct {
  int n;                        /* number of states */
  double *a;                    /* transition probability */
  double *pi;                   /* initial state */
  double (*emission) (void*, int, int); /* log emission probability, (observation, state_index) */
} dsmc_t;

double viterbi(dsmc_t *m, int t_max, void *o, int *q, int q_init, int verbose);

#endif /* _WZHMM_H */
