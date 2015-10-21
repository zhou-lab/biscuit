#include <stdlib.h>
#include <string.h>
#include "wzhmm.h"
#define _HMM_HUGE 10000000000.0
#define _HMM_SMALL 0.0000000001

/* input o, l. output alpha, scale, p */
double forward(dsmc_t *m, int t_max, void *o, double *alpha, double *scale) {

  int i, j;

  /* initialization */
  scale[0] = 0.0;
  for (i=0; i<m->n; ++i) {
    alpha[i] = m->pi[i]*m->emission(o, 0, i);
    scale[i] += alpha[i];
  }
  for (i=0; i<m->n; ++i) alpha[i] /= scale[i]; /* scale alpha */

  /* induction */
  int t;
  for (t=1; t<t_max; ++t) {
    scale[t] = 0.0;
    for (i=0; i<m->n; ++i) {    /* current state */
      double sum = 0.0;
      for (j=0; j<m->n; ++j)    /* previous state */
        sum += alpha[(t-1)*m->n+j]*m->a[j*m->n+i];

      alpha[t*m->n+i] = sum*m->emission(o, t, i);
      scale[t] += alpha[t*m->n+i];
    }
    for (i=0; i<m->n; ++i) alpha[t*m->n+i] /= scale[t]; /* scale alpha */
  }

  /* termination */
  double p=0;
  for (i=0; i<m->n; ++i) p+=alpha[t*m->n+i];
  return p;
}

/* comput beta */
void backward(dsmc_t *m, int t_max, void *o, double *beta, double *scale) {
  int i, j;
  int t;

  /* initialization */
  for (i=0; i<m->n; ++i) beta[(t_max-1)*m->n+i] = 1.0/scale[t_max-1];
  /* induction */
  for (t=m->n-1; t; --t) {
    for (i=0; i<m->n; ++i) {
      double sum = 0.0;
      for (j=0; j<m->n; ++j) sum += m->a[i*m->n+j] * m->emission(o, t, j)*beta[t*m->n+j];
      beta[(t-1)*m->n+i] = sum/scale[t-1];
    }
  }
}

void posterior_decoding(dsmc_t *m, int t_max, double *alpha, double *beta, int *q) {
  int i, t;

  double max_gamma;
  for (t=0; t<t_max; ++t) {
    for (i=0; i<m->n; ++i) {    /* current state */
      double gamma = alpha[t*m->n+i]*beta[t*m->n+i];
      if (!i || gamma > max_gamma) {
        max_gamma = gamma;
        q[t] = i;
      }
    }
  }
}

/* q_init can be arbitrarily set to 0 */
double viterbi(dsmc_t *m, int t_max, void *o, int *q, int q_init, int verbose) {

  int i,j;
  double *delta = calloc(m->n, sizeof(double));
  double *delta0 = calloc(m->n, sizeof(double));
  int *psi = calloc(t_max*m->n, sizeof(int));

  /* take log of transition probability */
  double *log_a = calloc(m->n*m->n, sizeof(double));
  for (i=0; i<m->n; ++i)
    for (j=0; j<m->n; ++j)
      log_a[i*m->n+j] = log(m->a[i*m->n+j]);

  /* initialization */
  for (i=0; i<m->n; ++i) {
    delta0[i] = log(m->pi[i]) + m->emission(o, 0, i);
    psi[i] = q_init;
  }

  /* induction */
  int t;
  for (t=1; t<t_max; ++t) {
    for (i=0; i<m->n; ++i) {    /* current state */
      double maxval=0.0; int maxvalind=0;
      for (j=0; j<m->n; ++j) {  /* previous state */
        double val = delta0[j] + log_a[j*m->n+i];
        if (j==0 || val > maxval) {
          maxval = val;
          maxvalind = j;
        }
      }
      if (verbose>6) {
        fprintf(stdout, "maxval: %1.3f\t%d\n", maxval, maxvalind);
      }
      delta[i] = maxval + m->emission(o, t, i);
      psi[t*m->n+i] = maxvalind;
    }
    memcpy(delta0, delta, m->n*sizeof(double));
  }

  /* termination
   * pick the largest delta as probability of the
   * the best parse and its index the state of the
   * last item */
  double p = 0.0;
  for (i=0; i<m->n; ++i) {
    if (i==0 || delta[i] > p) {
      p = delta[i];
      q[t_max-1] = i;
    }
  }

  /* state backtracking */
  for (t=t_max-1; t; --t)
    q[t-1] = psi[t*m->n+q[t]];

  if (verbose>1) {
    fprintf(stdout, "q (psi):");
    for (t=0; t<t_max; ++t)
      fprintf(stdout, "\t%d(%d)", q[t], psi[t]);
    fprintf(stdout, "\n");
  }

  free(log_a);
  free(delta); free(delta0); free(psi);
  return p;
}

void compute_gamma(dsmc_t *m, int t_max, double *alpha, double *beta, double *gamma) {
  int i;
  int t;
  for (t=0; t<t_max; ++t) {
    double denom = 0.0;
    for (i=0; i<m->n; ++i) {
      gamma[t*m->n+i] = alpha[t*m->n+i]*beta[t*m->n+i];
      denom += gamma[t*m->n+i];
    }
    for (i=0; i<m->n; ++i) 
      gamma[t*m->n+i] = gamma[t*m->n+i]/denom;
  }
}

void compute_xi(dsmc_t *m, int t_max, void *o, double *alpha, double *beta, double *xi) {
  int i,j;
  int t;
  for (t=0; t<t_max-1; ++t) {
    double sum=0.0;
    for (i=0; i<m->n; ++i) {
      for (j=0; j<m->n; ++j) {
        xi[t*m->n*m->n+i*m->n+j] = alpha[t*m->n+i]*beta[(t+1)*m->n+j]*m->a[i*m->n+j]*m->emission(o, t+1, j);
        sum += xi[t*m->n*m->n+i*m->n+j];
      }
    }

    for (i=0; i<m->n; ++i)
      for (j=0; j<m->n; ++j)
        xi[t*m->n*m->n+i*m->n+j] /= sum;
  }
}

void baum_welch(dsmc_t *m, int t_max, void *o) {

  int i,j;
  double *xi=calloc((t_max-1)*m->n*m->n, sizeof(double));
  double *gamma=calloc(t_max*m->n, sizeof(double));
  double *alpha = calloc(t_max*m->n, sizeof(double));
  double *beta = calloc(t_max*m->n, sizeof(double));
  double *scale = calloc(t_max, sizeof(double));
  double p, pp=0.0; double dp=0.001;
  int max_iter=10;
  for (i=0; i<max_iter; ++i) {
    p = forward(m, t_max, o, alpha, scale);
    backward(m, t_max, o, scale, &p);
    compute_gamma(m, t_max, alpha, beta, gamma);
    compute_xi(m, t_max, o, alpha, beta, xi);

    /* update pi */
    for (i=0; i<m->n; ++i) {
      m->pi[i] = gamma[i];
    }

    /* update a */
    for (i=0; i<m->n; ++i) {
      double sum_gamma = 0.0;
      int t;
      for (t=0; t<t_max; ++t) sum_gamma += gamma[m->n*t+i];

      for (j=0; j<m->n; ++j) {
        double sum_xi = 0.0;
        for (t=0; t<t_max; ++t) sum_xi += xi[t*m->n*m->n+i*m->n+j];
        m->a[i*m->n+j] = sum_gamma / sum_xi;
      }
    }

    if (i!=0 && p-pp < dp) break;
    pp = p;
  }
}


void naturalize(dsmc_t *m) {
  int i,j;
  /* make sure prior probability is higher than zero */
  for (i=0; i<m->n; ++i)
    if (m->pi[i]<=0) m->pi[i] = _HMM_SMALL;

  /* make sure transition probability is higher than zero */
  for (i=0; i<m->n; ++i) {
    for (j=0; j<m->n; ++j) {
      if (m->a[i*m->n+j]<=0) {
        m->a[i*m->n+j] = _HMM_SMALL;
      }
    }
  }
}

