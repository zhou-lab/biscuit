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
  double *a;                    /* transition probability, dim: n x n */
  double *pi;                   /* initial state */
  double (*emission) (void*, int, int); /* log emission probability, (observation, state_index) */
} dsmc_t;

static inline dsmc_t *init_dsmc(int n, double (*emission) (void*, int, int))  {
  dsmc_t *m = (dsmc_t*) calloc(1,sizeof(dsmc_t));
  m->n = n;
  m->a = calloc(m->n*m->n, sizeof(double));
  m->pi = calloc(m->n, sizeof(double));
  m->emission = emission;

  /* uniform prior probability */
  int i;
  for (i=0; i<m->n; ++i)
    m->pi[i] = 1.0/(double)m->n;

  return m;
}

static inline void free_dsmc(dsmc_t *m) {
  free(m->a); free(m->pi); free(m);
}

/* likelihood calculation - forward
 * decoding - viterbi, posterior
 * parameterization - baum-welch */

double viterbi(int *q, dsmc_t *m, int t_max, void *o, int q_init, int verbose);


#endif /* _WZHMM_H */
