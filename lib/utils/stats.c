#include "stats.h"

double beta_binomial(int k, int n, double a, double b) {

  double p = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) + lgamma(k+a) + lgamma(n-k+b) - lgamma(n+a+b) + lgamma(a+b) - lgamma(a) - lgamma(b);

  return p;
}

static double max( double a, double b ) { return a > b ? a : b; }

/* sum of numbers inside ln function,
 * i.e., ln_sum(ln(a), ln(b)) = ln(a+b) */
double ln_sum2(double ln1, double ln2){
  double I = max(ln1, ln2);
  return log(exp(ln1-I) + exp(ln2-I)) + I;
}

/* sum of numbers inside ln function, i.e.,
 * ln_sum(ln(a), ln(b), ln(c)) = ln(a+b+c) */
double ln_sum3(double ln1, double ln2, double ln3){
  double I = max(max(ln1, ln2), ln3);
  return log(exp(ln1-I) + exp(ln2-I) + exp(ln3-I)) + I;
}

/* sum of numbers inside ln function, i.e.,
 * ln_sum(ln(a), ln(b), ln(c), ln(d)) = ln(a+b+c+d) */
double ln_sum4(double ln1, double ln2, double ln3, double ln4){
  double I = max(max(max(ln1, ln2), ln3), ln4);
  return log(exp(ln1-I) + exp(ln2-I) + exp(ln3-I) + exp(ln4-I)) + I;
}

/* substract two numbers inside ln function, i.e.,
 * ln_sum(ln(a) - ln(b)) = ln(a-b) */
double ln_substract(double ln1, double ln2){
  double I = max(ln1, ln2);
  return log(exp(ln1-I) - exp(ln2-I)) + I;
}


static double binom_coeff(int m, int n, double p, double q) {
  double logcoeff = lgamma(m + n + 1.0);
  logcoeff -=  lgamma(n + 1.0) + lgamma(m + 1.0);
  logcoeff += m * log(p) + n * log(q);
  return exp(logcoeff);
}

/* survival p-value binomial distribution
 * s is the number of success,
 * n is the number of trial and
 * p is the error rate */

double binom_pval(int s, int n, double p) {
  double cdf = 1.0;
  register int i;
  for (i=0; i<s; i++) {
    cdf -= binom_coeff(i, n-i, p, 1-p);
  }
  return cdf;
}

/* Evaluate continued fraction for incomplete beta function
 * by modified Lentz's method
 *
 * copied from "Numerical Recipe in C",
 * 1992, Cambridge University Press */

static double beta_cf(double a, double b, double x) {
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;

  qab = a+b;
  qap = a+1.0;
  qam = a-1.0;
  c = 1.0;
  d = 1.0-qab * x / qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  h = d;
  for (m=1; m<=MAXIT; m++) {
    m2 = 2 * m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d = 1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c = 1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d = 1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c = 1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d = 1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) {
    perror("a or b too big, or MAXIT too small in betacf");
    abort();
  }
  return h;
}

static double beta_inc(double a, double b, double x) {
  double bt;
  if (x < 0.0 || x > 1.0) {
    fprintf(stderr, "a=%1.4g\tb=%1.4g\tx=%1.4g\n", a,b,x);
    perror("Bad x in routine betai");
    abort();
  }
  if (x == 0.0 || x == 1.0) bt=0.0;
  else /* Factors in front of the continued fraction. */
    bt = exp(lgamma(a+b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0)) /* Use continued fraction directly. */
    return bt * beta_cf(a,b,x) / a;
  else /* Use continued fraction after making the symmetry transformation. */
    return 1.0 - bt * beta_cf(b, a, 1.0-x) / b;
    
}



static double ln_beta_inc(double a, double b, double x) {
  double bt;
  if (x <= 0.0 || x >= 1.0) {
    fprintf(stderr, "a=%1.4g\tb=%1.4g\tx=%1.4g\n", a,b,x);
    perror("Bad x in ln_betai");
    abort();
  }

  /* Factors in front of the continued fraction. */
  bt = a * log(x) + b * log(1.0-x);
  
  if (x < (a+1.0)/(a+b+2.0)) /* Use continued fraction directly. */
    return bt + log(beta_cf(a,b,x)) - log(a);
  else  /* Use continued fraction after making the symmetry transformation. */
    /* printf("f: %1.4g\t%1.4g\t%1.4g\n", lgamma(a) + lgamma(b) - lgamma(a+b), bt, beta_cf(b, a, 1.0-x) / b); */
    return ln_substract(lgamma(a) + lgamma(b) - lgamma(a+b), bt + log(beta_cf(b, a, 1.0-x) / b));
}

#define pv(f) ((double)(f) * (1 - error) + (1 - (double)(f)) * error)

/* double */
/* pval2qual(double pval) { */
/*   return - 10 * log10(max(pval, pow(10, -25.5))); */
/* } */

double
pval2qual(double pval) {
  int qual = (int) (- 10 * log10(max(pval, pow(10, -30.0))) + 0.499);
  if (qual > 255) return 255;
  else return qual;
}

double
varcall_pval(int kr, int kv, double error, double mu, double contam) {

  double U, V;
  if (contam == 0.0) {
    U = pow(pv(0), kv) * pow(1-pv(0), kr) * (1-mu);
    V = (beta_inc(kv+1, kr+1, pv(1.0)) - beta_inc(kv+1, kr+1, pv(0))) * mu;
  } else if (contam >= 0.0) {
    U = (beta_inc(kv+1, kr+1, pv(contam)) - beta_inc(kv+1, kr+1, pv(0))) * (1-mu);
    V = (beta_inc(kv+1, kr+1, pv(1.0)) - beta_inc(kv+1, kr+1, pv(0))) * mu * contam;
  } else {
    perror("Contamination extent cannot be negative.");
    abort();
  }
  
  return U/(U+V);
}

double
ref_lnlik(int kr, int kv, double error, double contam) {
  if (contam == 0.0) {
    return ln_binom_kernel(pv(0.0), kv, kr) - log(1-2*error) + lgamma(kr+kv+1) - lgamma(kv+1) - lgamma(kr+1);
  } else if (contam > 0.0) {
    return ln_beta_incdiff_kernel(pv(0.0), pv(contam), kv+1, kr+1) - log(1-2*error) + lgamma(kr+kv+1) - lgamma(kv+1) - lgamma(kr+1) - log(contam);
  } else {
    perror("Contamination extent cannot be negative.");
    abort();
  }
}

double
alt_lnlik(int kr, int kv, double error) {
  return ln_beta_incdiff_kernel(pv(0.0), pv(1.0), kv+1, kr+1) - log(1-2*error) + lgamma(kr+kv+1) - lgamma(kv+1) - lgamma(kr+1);
}


double
somatic_posterior(int kr_t, int kv_t, int kr_n, int kv_n, double error, double mu, double mu_somatic, double contam) {

  double prob_m00 = ref_lnlik(kr_n, kv_n, error, contam) + \
    ref_lnlik(kr_t, kv_t, error, contam);
  double prob_m01 = ref_lnlik(kr_n, kv_n, error, contam) + \
    alt_lnlik(kr_t, kv_t, error) + log(mu_somatic);
  double prob_m10 = alt_lnlik(kr_n, kv_n, error) +  \
    ref_lnlik(kr_t, kv_t, error, contam) + log(mu_somatic);
  double prob_m11 = alt_lnlik(kr_n, kv_n, error) +  \
    alt_lnlik(kr_t, kv_t, error) + log(mu);
  double prob_d = ln_sum4(prob_m00, prob_m01, prob_m10, prob_m11);

#ifdef DEBUGSTATS

  printf("t %d, %d, ref_lnlik: %1.7f\n", kr_t, kv_t, ref_lnlik(kr_t, kv_t, error, contam));
  printf("t %d, %d, alt_lnlik: %1.7f\n", kr_t, kv_t, alt_lnlik(kr_t, kv_t, error));
  printf("n %d, %d, ref_lnlik: %1.7f\n", kr_n, kv_n, ref_lnlik(kr_n, kv_n, error, contam));
  printf("n %d, %d, alt_lnlik: %1.7f\n", kr_n, kv_n, alt_lnlik(kr_n, kv_n, error));

  printf("m00: %1.7f\n", prob_m00);
  printf("m01: %1.7f\n", prob_m01);
  printf("m10: %1.7f\n", prob_m10);
  printf("m11: %1.7f\n", prob_m11);
  printf("d: %1.7f\n", prob_d);

#endif // DEBUGSTATS

  return 1-exp(prob_m01 - prob_d);
}

double
genotype_prior_HWE(Genotype genotype, double allele_freq) {
  if (genotype == HOMOVAR) {
    return allele_freq * allele_freq;
  } else if (genotype == HET) {
    return 2 * allele_freq * (1-allele_freq);
  } else {
    return 1.0 - 2 * allele_freq * (1-allele_freq) - allele_freq * allele_freq;
  }
}

double
ln_binom_kernel(double p, int a, int b) {
  /* prevent numerical underflow */
  if (p < FPMIN)
    p = FPMIN;
  if (p > 1-FPMIN)
    p = 1-FPMIN;
  return log(p) * a + log(1-p) * b;
}

double
ln_beta_incdiff_kernel(double p1, double p2, int a, int b) {
  /* note that p2 must be greater than p1 */
  /* symmetric transform if both p1 and p2 are on the large side */
  if (p1 > (a+1.0)/(a+b+2.0) && p2 > (a+1.0)/(a+b+2.0)) {
    double tmp_dbl;
    tmp_dbl = p1;
    p1 = 1 - p2;
    p2 = 1 - tmp_dbl;
    int tmp_int;
    tmp_int = a;
    a = b;
    b = tmp_int;
  }
  return log(1 - exp(ln_beta_inc(a,b,p1) - ln_beta_inc(a,b,p2))) + ln_beta_inc(a,b,p2);
}

double
somatic_lnlik(int kr, int kv, double error) {
  return ln_beta_incdiff_kernel(pv(0.0), pv(1.0), kv+1, kr+1) - log(1-2*error) + lgamma(kv+kr+1) - lgamma(kv+1) - lgamma(kr+1);
}

double
inconsist_score(int kr_tumor, int kv_tumor,
		int kr_normal, int kv_normal,
		double mu, double error) {
  int kv = kv_normal + kv_tumor;
  int kr = kr_normal + kr_tumor;

  double coeff1 = lgamma(kv_tumor+kr_tumor+1) - lgamma(kv_tumor+1) - lgamma(kr_tumor+1);
  double coeff2 = lgamma(kv_normal+kr_normal+1) - lgamma(kv_normal+1) - lgamma(kr_normal+1);
  double consist_lnlik = ln_beta_incdiff_kernel(pv(0.0), pv(1.0), kv+1, kr+1) - log(1-2*error) + coeff1 + coeff2;
  double inconsist_lnlik = ln_sum2(somatic_lnlik(kr_tumor, kv_tumor, error), somatic_lnlik(kr_normal, kv_normal, error));
  
  return - consist_lnlik - log(1-mu) + ln_sum2(consist_lnlik + log(1-mu), inconsist_lnlik + log(mu));
}

/* degeneracy: reference (dr=1), variant (dv=0) */
double
genotype_lnlik(Genotype genotype, int kr, int kv, double error, double contam) {
  double lnlik;
  switch (genotype) {
  case HOMOREF:
    if (contam == 0.0) {
      lnlik = ln_binom_kernel(pv(0.0), kv, kr);
    } else if (contam >= 0.0) {
      lnlik = ln_beta_incdiff_kernel(pv(0.0), pv(contam), kv+1, kr+1) - log(contam) - log(1-2*error);
    } else {
      perror("Contamination extent cannot be negative.");
      abort();
    }
    break;

  case HET:
    if (contam == 0.0) {
      lnlik = ln_binom_kernel(pv(0.5), kv, kr);
    } else if (contam >= 0.0) {
      lnlik = ln_beta_incdiff_kernel(pv(0.5-contam), pv(0.5+contam), kv+1, kr+1) - log(2*contam) - log(1-2*error);
    } else {
      perror("Contamination extent cannot be negative.");
      abort();
    }
    break;

  case HOMOVAR:
    if (contam == 0.0) {
      lnlik = ln_binom_kernel(pv(1.0), kv, kr);
    } else if (contam >= 0.0) {
      lnlik = ln_beta_incdiff_kernel(pv(1.0-contam), pv(1.0), kv+1, kr+1) - log(contam)- log(1-2*error);
    } else {
      perror("Contamination extent cannot be negative.");
      abort();
    }
    break;

  default:
    perror("Genotype not recognized\n");
    abort();
  }

  return lnlik + lgamma(kv+kr+1) - lgamma(kv+1) - lgamma(kr+1);
}

/**** Fisher's exact test ****/
/* log\binom{n}{k} */
static double lbinom(int n, int k)
{
  if (k == 0 || n == k) return 0;
  return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

/* n11  n12  | n1_
 * n21  n22  | n2_
 *-----------+----
 * n_1  n_2  | n
 */

/* hypergeometric distribution */
static double hypergeo(int n11, int n1_, int n_1, int n)
{
  return exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
}

typedef struct {
  int n11, n1_, n_1, n;
  double p;
} hgacc_t;

// incremental version of hypergenometric distribution
static double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux) {
  if (n1_ || n_1 || n) {
    aux->n11 = n11; aux->n1_ = n1_; aux->n_1 = n_1; aux->n = n;
  } else { // then only n11 changed; the rest fixed
    if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
      if (n11 == aux->n11 + 1) { // incremental
        aux->p *= (double)(aux->n1_ - aux->n11) / n11
          * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
        aux->n11 = n11;
        return aux->p;
      }
      if (n11 == aux->n11 - 1) { // incremental
        aux->p *= (double)aux->n11 / (aux->n1_ - n11)
          * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
        aux->n11 = n11;
        return aux->p;
      }
    }
    aux->n11 = n11;
  }
  aux->p = hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
  return aux->p;
}

double fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two) {
  int i, j, max, min;
  double p, q, left, right;
  hgacc_t aux;
  int n1_, n_1, n;

  n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
  max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
  min = n1_ + n_1 - n;    // not sure why n11-n22 is used instead of min(n_1,n1_)
  if (min < 0) min = 0; // min n11, for left tail
  *two = *_left = *_right = 1.;
  if (min == max) return 1.; // no need to do test
  q = hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
  // left tail
  p = hypergeo_acc(min, 0, 0, 0, &aux);
  for (left = 0., i = min + 1; p < 0.99999999 * q && i<=max; ++i) // loop until underflow
    left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
  --i;
  if (p < 1.00000001 * q) left += p;
  else --i;
  // right tail
  p = hypergeo_acc(max, 0, 0, 0, &aux);
  for (right = 0., j = max - 1; p < 0.99999999 * q && j>=0; --j) // loop until underflow
    right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
  ++j;
  if (p < 1.00000001 * q) right += p;
  else ++j;
  // two-tail
  *two = left + right;
  if (*two > 1.) *two = 1.;
  // adjust left and right
  if (abs(i - n11) < abs(j - n11)) right = 1. - left + q;
  else left = 1.0 - right + q;
  *_left = left; *_right = right;
  return q;
}

#ifdef STAT_TESTMAIN
int main(int argc, char **argv) {
  /* printf("posterior: %1.7f\n", somatic_posterior(635, 633, 373, 1, 0.001, 0.001, 0.001, 0.001)); */
  /* printf("%1.7f\n", pval2qual(somatic_posterior(635, 20, 373, 0, 0.001, 0.001, 0.001, 0.001))); */
  printf("%1.7f\n", pval2qual(varcall_pval(827, 5, 0.001, 0.001, 0.01)));
  return 0;
}
#endif
