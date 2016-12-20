typedef struct {
  int verbose;
  /* for emission probability */
  int a;
  int b;
  /* for transition probability */
  double p;
  double q;
} conf_t;

int main_asm(int argc, char *argv[]) {

  conf_t conf = {.verbose=6, .a=10, .b=10, .p=0.5, .q=0.5};

  methbed_t *in = calloc(1, sizeof(methbed_t));
  int c, i; int collapse=0; char *out_fn=0;
  while ((c = getopt(argc, argv, "i:")) >= 0) {
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
      fprintf(stderr, "Usage: biscuit asm [options] pairwise_epiread \n");
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

  if (optind >= argc) {
    fprintf(stderr, "Please provide input epiread.\n");
    fflush(stderr);
  }

  char *input_fn = argv[optind];
  read_tsv(input_fn);

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
