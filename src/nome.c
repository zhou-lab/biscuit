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
  char *vcf_fn;
  char *bed_fn;
  gzFile FH;
  meth_obs1_t next;             /* next object */
  kstring_t nextchrom;          /* next chromosome */
  int8_t isvcf;
} pileup_t;

double meth_emission(void *d, int t, int state_index) {
  meth_obs1_t *ob1 = (meth_obs1_t*) d;
  if (state_index == 1) 
    return beta_binomial(ob1[t].ret, ob1[t].cov, 1, 10);
  else
    return beta_binomial(ob1[t].ret, ob1[t].cov, 10, 1);
}

static int vcf_parse1(char *line, meth_obs1_t *ob, char *chrom) {

  char *tok;
  char *linerest=0, *fieldrest=0;

  tok=strtok_r(line, "\t", &linerest);
  strcpy(chrom, tok);

  tok=strtok_r(NULL, "\t", &linerest);
  ensure_number(tok);
  int64_t pos = atoi(tok);
  
  /* ID */
  tok=strtok_r(NULL, "\t", &linerest);

  /* REF */
  tok=strtok_r(NULL, "\t", &linerest);

  /* ALT */
  tok=strtok_r(NULL, "\t", &linerest);

  /* QUAL */
  tok=strtok_r(NULL, "\t", &linerest);

  /* PASS FILTER */
  tok=strtok_r(NULL, "\t", &linerest);

  /* INFO */
  tok=strtok_r(NULL, "\t", &linerest);
  int i; int is_gch=0;
  for (i=0; i<(signed)strlen(tok)-2; ++i) {
    if (strncmp(tok+i,"N5=",3)==0) {
      if (tok[i+4]=='G' && tok[i+5]=='C' && tok[i+6]!='G') {
        is_gch=1;
        break;
      }
    }
  }

  if (!is_gch) return 0;

  /* FORMAT */
  tok=strtok_r(NULL, "\t", &linerest);
  int cv_index=-1, bt_index=-1;
  char *field;
  field=strtok_r(tok, ":", &fieldrest);
  i=0;
  while(field) {
    if (field[0]=='C' && field[1]=='V') {
      cv_index = i;
    }
    if (field[0]=='B' && field[1]=='T') {
      bt_index = i;
    }
    ++i;
    field=strtok_r(NULL, ":", &fieldrest);
  }

  if (cv_index<0 || bt_index<0) return 0;

  /* FORMAT content */
  tok = strtok_r(NULL, "\t", &linerest);
  i=0;
  int coverage=-1; double beta=-1.0;
  field=strtok_r(tok, ":", &fieldrest);
  while (field) {
    if (i==cv_index) {
      ensure_number(field);
      coverage=atoi(field);
    }
    if (i==bt_index) {
      ensure_number(field);
      beta=atof(field);
    }
    ++i;
    field=strtok_r(NULL, ":", &fieldrest);
  }

  if (coverage<0||beta<0) return 0;

  ob->pos = pos;
  ob->cov = coverage;
  ob->ret = (int) (coverage*beta);

  return 1;
}

meth_obs1_v *pileup_get_chrom1(pileup_t *in, char **chrom) {

  meth_obs1_v *obs=init_meth_obs1_v(10000);
  meth_obs1_t *ob = try_next_meth_obs1_v(obs);

  *chrom = 0;
  if (in->nextchrom.l) {
    *ob = in->next;
    commit_next_meth_obs1_v(obs);
    ob = try_next_meth_obs1_v(obs);
    *chrom = strdup(in->nextchrom.s); /* record current chrom once */
  }

  char ch[1000];
  kstring_t str;
  str.s = 0; str.l = str.m = 0;

  int i=0;
  while (1) {
    int c=gzgetc(in->FH);
    if (c=='\n' || c==EOF) {
      if (*chrom && i%100000==0) {
        fprintf(stderr, "\r%s\t%d\t%zu", *chrom, i, obs->size);
        fflush(stderr);
      }
      ++i;
      if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) {
        if (vcf_parse1(str.s, ob, ch)) {
          if (!in->nextchrom.l) {
            kputs(ch, &in->nextchrom);
            *chrom = strdup(in->nextchrom.s); /* record current chrom once */
          } else if (strcmp(ch, in->nextchrom.s) != 0) {
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
      if (c==EOF) break;
    } else {
      kputc(c, &str);
    }
  }

  return obs;
}

void pileup_open(pileup_t *in) {

  memset(&in->next, 0, sizeof(meth_obs1_t));
  in->nextchrom.s = 0;
  in->nextchrom.l = in->nextchrom.m = 0;
  if (in->vcf_fn) {
    in->FH = gzopen(in->vcf_fn,"r");
    in->isvcf = 1;
  } else if (in->bed_fn) {
    in->FH = gzopen(in->bed_fn,"r");
    in->isvcf = 0;
  } else {
    fprintf(stderr, "[%s:%d] Please provide pileup input.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
}

void pileup_close(pileup_t *in) {
  free(in->nextchrom.s);
  gzclose(in->FH);
}

int main_nome(int argc, char *argv[]) {

  conf_t conf = {.verbose=6};

  pileup_t in = {.vcf_fn=0, .FH=0, .bed_fn=0};
  int c, i;
  while ((c = getopt(argc, argv, "V:i:h")) >= 0) {
    switch (c) {
    case 'V': conf.verbose = atoi(optarg); break;
    case 'i': in.vcf_fn = optarg; break;
    case 'b': in.bed_fn = optarg; break;
    case 'h': {
      fprintf(stderr, "\n");
      fprintf(stderr, "Usage: biscuit nome [options] -i in.vcf \n");
      fprintf(stderr, "Input options:\n");
      fprintf(stderr, "     -i FILE   input vcf, coordinates-sorted\n");
      fprintf(stderr, "     -b FILE   bed file of GpC methylation coordinate-sorted\n");
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

  pileup_open(&in);
  char *chrom;
  meth_obs1_v *obs=pileup_get_chrom1(&in, &chrom);

  if (conf.verbose>3) {
    for (i=0; i<(signed)obs->size; ++i) {
      meth_obs1_t *o = ref_meth_obs1_v(obs,i);
      fprintf(stdout, "%s\t%"PRId64"\t%d\t%d\n", chrom, o->pos, o->cov, o->ret);
    }
  }
  
  pileup_close(&in);

  /* make a 2-state hmm */
  dsmc_t *m = (dsmc_t*) calloc(1,sizeof(dsmc_t));
  m->n = 2;
  m->a = calloc(m->n*m->n, sizeof(double));
  m->pi = calloc(m->n, sizeof(double));
  m->emission = meth_emission;

  m->a[0*2] = 0.5;
  m->a[0*2+1] = 0.5;
  m->a[1*2] = 0.5;
  m->a[1*2+1] = 0.5;

  m->pi[0] = 0.5;
  m->pi[1] = 0.5;
  
  /* meth_obs1_t obs[1000]; */
  for (i=0; i<100; ++i) {
    if ((i/10)%2==1) {
      obs->buffer[i].cov = 30;
      obs->buffer[i].ret = 0;
    } else {
      obs->buffer[i].cov = 10;
      obs->buffer[i].ret = 8;
    }
  }

  /* double p; */
  /* double *alpha = calloc(100, sizeof(double)); */
  /* double *scale = calloc(100, sizeof(double)); */
  /* p=forward(m, 100, obs, alpha, scale); */
  int *q = calloc(100, sizeof(int));
  viterbi(m, 100, obs->buffer, q, 0, conf.verbose);

  free(q);
  free(m->a); free(m->pi); free(m);
  
  return 0;
}
