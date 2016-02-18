#include "wzhmm.h"
#include "stats.h"
#include "kstring.h"
#include "wvec.h"
#include "wstr.h"

typedef struct {
  int verbose;
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

double meth_emission(void *d, int t, int state_index) {
  meth_obs1_t *ob1 = (meth_obs1_t*) d;
  if (state_index == 1) 
    return beta_binomial(ob1[t].ret, ob1[t].cov, 1, 10);
  else
    return beta_binomial(ob1[t].ret, ob1[t].cov, 10, 1);
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
        fprintf(stderr, "\r%s\t%d\t%zu", in->chrm, i, obs->size);
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

  conf_t conf = {.verbose=6};

  methbed_t *in = calloc(1, sizeof(methbed_t));
  int c, i;
  while ((c = getopt(argc, argv, "V:b:h")) >= 0) {
    switch (c) {
    case 'V': conf.verbose = atoi(optarg); break;
    case 'b': in->bed = optarg; break;
    case 'h': {
      fprintf(stderr, "\n");
      fprintf(stderr, "Usage: biscuit ndr [options] -b in.bed \n");
      fprintf(stderr, "Input options:\n");
      fprintf(stderr, "     -b FILE   bed file of GpC retention, coordinate-sorted\n");
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

  methbed_open(in);

  while (1) {

    meth_obs1_v *obs=methbed_get_chrom1(in);
    if (!obs->size) break;

    if (conf.verbose>3) {
      for (i=0; i<(signed)obs->size; ++i) {
        meth_obs1_t *o = ref_meth_obs1_v(obs,i);
        fprintf(stdout, "%s\t%"PRId64"\t%d\t%d\n", in->chrm, o->pos, o->cov, o->ret);
      }
    }
  
    /* make a 2-state hmm */
    dsmc_t *m = init_dsmc(2, meth_emission);

    m->a[0*2] = 0.5;
    m->a[0*2+1] = 0.5;
    m->a[1*2] = 0.5;
    m->a[1*2+1] = 0.5;

    /* /\* fake an observation vector, meth_obs1_t obs[1000]; *\/ */
    /* for (i=0; i<100; ++i) { */
    /*   if ((i/10)%2==1) { */
    /*     obs->buffer[i].cov = 30; */
    /*     obs->buffer[i].ret = 0; */
    /*   } else { */
    /*     obs->buffer[i].cov = 10; */
    /*     obs->buffer[i].ret = 8; */
    /*   } */
    /* } */

    int *q = calloc(obs->size, sizeof(int));
    viterbi(q, m, obs->size, obs->buffer, 0, conf.verbose);

    unsigned j;
    for (j=0; j<obs->size; ++j) {
      meth_obs1_t *o = ref_meth_obs1_v(obs,j);
      fprintf(stdout, "%s\t%"PRId64"\t%d\n", in->chrm, o->pos, q[j]);
      fflush(stdout);
    }

    free_meth_obs1_v(obs);
    free(q);
    free_dsmc(m);
  }

  methbed_close(in);
  free_methbed(in);

  return 0;
}
