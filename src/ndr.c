/**
 * Find nucleosome-depleted region.
 *
 * The MIT License (MIT)
 * Copyright (c) 2016 Wanding.Zhou@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
**/

#include "wzhmm.h"
#include "stats.h"
#include "kstring.h"
#include "wvec.h"
#include "wzmisc.h"

typedef struct {
  int verbose;
  /* for emission probability */
  int a;
  int b;
  /* for transition probability */
  double p;
  double q;
} conf_t;

/* observation of methylation */
typedef struct {
  int64_t pos;
  int cov;                      /* coverage */
  int ret;                      /* retention */
} meth_obs1_t;

DEFINE_VECTOR(meth_obs1_v, meth_obs1_t);

typedef struct {
  char *bed;
  gzFile FH;
  char *chrm;
  meth_obs1_t next;             /* next object */
  kstring_t nextchrom;          /* chromosome of next object */
} methbed_t;

void methbed_close(methbed_t *m) {
  gzclose(m->FH);
}

void free_methbed(methbed_t *m) {
  free(m->chrm);
  free(m->nextchrom.s);
  free(m);
}

double meth_emission(void *d, int t, int state_index, void *_c) {
  meth_obs1_t *ob1 = (meth_obs1_t*) d;
  conf_t *c = (conf_t*) _c;
  if (state_index == 1)
    return beta_binomial(ob1[t].ret, ob1[t].cov, 1, c->b);
  else
    return beta_binomial(ob1[t].ret, ob1[t].cov, c->a, 1);
}

int methbed_parse1(char *line, meth_obs1_t *ob, char *chrom) {

  char *tok;

  tok=strtok(line, "\t");
  strcpy(chrom, tok);

  /* start */
  tok=strtok(NULL, "\t");

  /* end */
  tok=strtok(NULL, "\t");
  ensure_number(tok);
  ob->pos = atoi(tok);          /* 1-based */

  /* beta */
  tok=strtok(NULL, "\t");
  ensure_number(tok);
  double beta = atof(tok);

  /* coverage */
  tok=strtok(NULL, "\t");
  ensure_number(tok);
  ob->cov = atoi(tok);

  ob->ret = (int) (ob->cov*beta);

  return 1;
}

meth_obs1_v *methbed_get_chrom1(methbed_t *in) {

  meth_obs1_v *obs=init_meth_obs1_v(10000);
  meth_obs1_t *ob = try_next_meth_obs1_v(obs);

  in->chrm = 0;
  if (in->nextchrom.l) {        /* get next object from last run */
    *ob = in->next;
    commit_next_meth_obs1_v(obs);
    ob = try_next_meth_obs1_v(obs);
    in->chrm = strdup(in->nextchrom.s); /* record current chrom once */
  }

  char ch[1000];
  kstring_t str;
  str.s = 0; str.l = str.m = 0;

  int i=0;
  while (1) {
    int c=gzgetc(in->FH);
    if (c=='\n' || c==EOF) {

      if (in->chrm && i%100000==0) {
        fprintf(stderr, "\r[%s:%d] %s\t%d\t%zu\033[K", __func__, __LINE__, in->chrm, i, obs->size);
        fflush(stderr);
      }
      ++i;

      if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=3) {
        if (methbed_parse1(str.s, ob, ch)) {
          if (!in->nextchrom.l) { /* first time read */
            kputs(ch, &in->nextchrom);
            in->chrm = strdup(in->nextchrom.s); /* record current chrom once */
          } else if (strcmp(ch, in->nextchrom.s) != 0) { /* next chromosome encountered */
            in->next = *ob;
            in->nextchrom.l = 0;
            kputs(ch, &in->nextchrom);
            break;
          }
          commit_next_meth_obs1_v(obs);
          ob = try_next_meth_obs1_v(obs);
        }
      }
      str.l = 0;                /* clean line */
      if (c==EOF) {
        memset(&in->next, 0, sizeof(meth_obs1_t));
        in->nextchrom.l = 0;
        break;
      }
    } else {
      kputc(c, &str);
    }
  }

  free(str.s);

  return obs;
}

void methbed_open(methbed_t *in) {

  memset(&in->next, 0, sizeof(meth_obs1_t));
  in->nextchrom.s = 0;
  in->nextchrom.l = in->nextchrom.m = 0;
  if (in->bed) {
    if (strcmp(in->bed, "-") == 0) {
      in->FH = gzdopen(fileno(stdin), "r");
    } else {
      in->FH = gzopen(in->bed,"r");
      if (!in->FH) {
        fprintf(stderr, "[%s:%d] Cannot open %s\n", __func__, __LINE__, in->bed);
        fflush(stderr);
        exit(1);
      }
    }
  } else {
    in->FH = gzdopen(fileno(stdin), "r");
  }
}

int main_ndr(int argc, char *argv[]) {

  conf_t conf = {.verbose=6, .a=10, .b=10, .p=0.5, .q=0.5};

  methbed_t *in = calloc(1, sizeof(methbed_t));
  int c, i; int collapse=0; char *out_fn=0;
  while ((c = getopt(argc, argv, "V:b:o:A:B:P:Q:ch")) >= 0) {
    switch (c) {
    case 'V': conf.verbose = atoi(optarg); break;
    case 'b': in->bed = optarg; break;
    case 'o': out_fn = optarg; break;
    case 'A': conf.a = atoi(optarg); break;
    case 'B': conf.b = atoi(optarg); break;
    case 'P': conf.p = atof(optarg); break;
    case 'Q': conf.q = atof(optarg); break;
    case 'c': collapse = 1; break;
    case 'h': {
      fprintf(stderr, "\n");
      fprintf(stderr, "Usage: biscuit ndr [options] -b in.bed \n");
      fprintf(stderr, "Input options:\n");
      fprintf(stderr, "     -b FILE   bed file of GpC retention, coordinate-sorted\n");
      fprintf(stderr, "     -c        collapse into regions\n");
      fprintf(stderr, "     -A INT    parameter a in beta binomial emission [%d]\n", conf.a);
      fprintf(stderr, "     -B INT    parameter b in beta binomial emission [%d]\n", conf.b);
      fprintf(stderr, "     -P FLOAT  parameter p in transition probability [%1.2f]\n", conf.p);
      fprintf(stderr, "     -Q FLOAT  parameter q in transition probability [%1.2f]\n", conf.q);
      fprintf(stderr, "     -o FILE   output file\n");
      fprintf(stderr, "     -V INT    verbose level [%d].\n", conf.verbose);
      fprintf(stderr, "     -h        this help.\n");
      fprintf(stderr, "\n");
      return 1;
    }
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      fflush(stderr);
      exit(1);
      break;
    }
  }

  if (!in->bed) {
    fprintf(stderr, "[%s:%d] Error, no bed input of GpC retention.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }

  FILE *out;
  if (out_fn)
    out = fopen(out_fn, "w");
  else
    out = stdout;

  methbed_open(in);

  while (1) {

    meth_obs1_v *obs=methbed_get_chrom1(in);
    if (!obs->size) {
      free_meth_obs1_v(obs);
      break;
    }

    if (conf.verbose>3) {
      for (i=0; i<(signed)obs->size; ++i) {
        meth_obs1_t *o = ref_meth_obs1_v(obs,i);
        fprintf(out, "%s\t%"PRId64"\t%d\t%d\n", in->chrm, o->pos, o->cov, o->ret);
      }
    }
  
    /* make a 2-state hmm */
    dsmc_t *m = init_dsmc(2, meth_emission, &conf);

    m->a[0*2] = conf.p;
    m->a[0*2+1] = 1.0 - conf.p;
    m->a[1*2] = conf.q;
    m->a[1*2+1] = 1.0 - conf.q;

    int *q = calloc(obs->size, sizeof(int));
    viterbi(q, m, obs->size, obs->buffer, 0, conf.verbose);

    unsigned j;
    if (collapse) {
      int64_t beg=-1, end;
      for (j=0; j<obs->size; ++j) {
        meth_obs1_t *o = ref_meth_obs1_v(obs,j);
        if (q[j] == 1) {
          end = o->pos;
          if (beg<0)            /* start a new region */
            beg = o->pos-1;
        } else {
          if (beg>0) {          /* end a region */
            fprintf(out, "%s\t%"PRId64"\t%"PRId64"\n", in->chrm, beg, end);
            beg = -1;
          }
        }
      }
      if (beg>=0)
        fprintf(out, "%s\t%"PRId64"\t%"PRId64"\n", in->chrm, beg, end);
      fflush(out);
    } else {
      for (j=0; j<obs->size; ++j) {
        meth_obs1_t *o = ref_meth_obs1_v(obs,j);
        fprintf(out, "%s\t%"PRId64"\t%"PRId64"\t%d\n", in->chrm, o->pos-1, o->pos, q[j]);
      }
      fflush(out);
    }


    free_meth_obs1_v(obs);
    free(q);
    free_dsmc(m);
    free(in->chrm);
  }
  fprintf(stderr, "\r[%s:%d] Done\033[K\n", __func__, __LINE__);
  fflush(stderr);

  methbed_close(in);
  free_methbed(in);

  if (out_fn) fclose(out);

  return 0;
}
