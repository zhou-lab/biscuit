#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "wstr.h"

typedef struct conf_t {
  char target[5];
  int mincov;
  int verbose;
} conf_t;

typedef struct bed1_t {
  char *chrm;
  int64_t pos;
  int64_t end;
  char ref;
  double beta;
  int cov;
} bed1_t;

bed1_t *init_bed1() {
  bed1_t *b=calloc(1, sizeof(bed1_t));
  b->chrm=0;
  return b;
}

void free_bed1(bed1_t *b) {
  free(b);
}

#define ET_C   0x1
#define ET_CG  0x2
#define ET_HCG 0x4
#define ET_GCH 0x8

static int vcf_parse1(char *line, bed1_t *b, uint8_t *et) {

  char *tok;
  char *linerest=0, *fieldrest=0;

  /* CHROM */
  tok=strtok_r(line, "\t", &linerest);
  b->chrm = realloc(b->chrm, strlen(tok)+1);
  strcpy(b->chrm, tok);

  /* POS */
  tok=strtok_r(NULL, "\t", &linerest);
  ensure_number(tok);
  b->pos = atoi(tok); b->end = atoi(tok);
  
  /* ID */
  tok=strtok_r(NULL, "\t", &linerest);

  /* REF */
  tok=strtok_r(NULL, "\t", &linerest);
  b->ref = tok[0];

  /* ALT */
  tok=strtok_r(NULL, "\t", &linerest);

  /* QUAL */
  tok=strtok_r(NULL, "\t", &linerest);

  /* PASS FILTER */
  tok=strtok_r(NULL, "\t", &linerest);

  /* INFO */
  tok=strtok_r(NULL, "\t", &linerest);
  *et=0;
  int i;
  for (i=0; i<(signed)strlen(tok)-2; ++i) {
    if (strncmp(tok+i,"N5=",3)==0) {
      if (tok[i+5]=='C') {
        *et |= ET_C;
        if (tok[i+4]=='G' && tok[i+6]!='G') *et |= ET_GCH;
        if (tok[i+6]=='G') {
          *et |= ET_CG;
          if (tok[i+4]!='G') *et |= ET_HCG;
        }
      }
      break;
    }
  }

  /* not any of the target */
  if (!*et) return 0;

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

  /* no coverage or beta info */
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

  /* no coverage or beta info */
  if (coverage<0||beta<0) return 0;

  b->cov = coverage;
  b->beta = beta;

  return 1;
}

void vcf2cg(gzFile FH, conf_t *conf) {

  bed1_t *b=init_bed1();        /* base */
  bed1_t *p=init_bed1();        /* previous base */

  kstring_t str;
  str.s = 0; str.l = str.m = 0;

  int et;
  int p_valid=0;
  while (1) {
    int c=gzgetc(FH);
    if (c=='\n' || c==EOF) {
      int merged=0;
      if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) {
        int parse_ok = vcf_parse1(str.s, b, &et);
        if (parse_ok && p_valid) {
          /* printf("%d\t%d\t%d\t%c\t%c\t%s\t%s\n", et, p->pos, b->pos, p->ref, b->ref, p->chrm, b->chrm); */
          if ((et&ET_CG) && p->pos+1 == b->pos && p->ref == 'C' && b->ref == 'G' && strcmp(p->chrm, b->chrm)==0) { /* merge CG */
            p->end++;
            /* fprintf(stdout, "merge%s\t%d\t%d\t%1.2f\t%d\n",b->chrm,b->pos,b->pos,b->beta,b->cov); */
            p->beta = (double)(p->beta*p->cov+b->beta*b->cov)/(p->cov+b->cov);
            p->cov += b->cov;
            merged=1;
          }
        }
        if (p_valid) {
          if (p->cov >= conf->mincov) {
            fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%1.3f\n", p->chrm, p->pos-1, p->end, p->beta);
          }
          p_valid = 0;
        }

        if (!merged && parse_ok && (et&ET_CG)) {
          /* fprintf(stdout, "%s\t%d\t%d\t%1.2f\t%d\n",b->chrm,b->pos,b->pos,b->beta,b->cov); */
          bed1_t *tmp;
          tmp = p; p = b; b = tmp;
          p_valid=1;
        }
      }
      str.l = 0;                /* clean line */
      if (c==EOF) break;
    } else {
      kputc(c, &str);
    }
  }
  if (p_valid) {
    if (p->cov >= conf->mincov) {
      fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%1.3f\n", p->chrm, p->pos-1, p->end, p->beta);
    }
  }
}

void vcf2hcg(gzFile FH, conf_t *conf) {

  bed1_t *b=init_bed1();
  bed1_t *p=init_bed1();

  kstring_t str;
  str.s = 0; str.l = str.m = 0;
  uint8_t et;
  int p_valid=0;
  while (1) {
    int c=gzgetc(FH);
    if (c=='\n' || c==EOF) {
      int merged=0;
      if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) {
        int parse_ok = vcf_parse1(str.s, b, &et);
        if (parse_ok && p_valid) {
          if ((et&ET_HCG) && p->pos+1 == b->pos && p->ref == 'C' && b->ref == 'G' && strcmp(p->chrm, b->chrm)==0) { /* merge CG */
            p->end++;
            p->beta = (double)(p->beta*p->cov+b->beta*b->cov)/(p->cov+b->cov);
            p->cov += b->cov;
            merged=1;
          }
        }
        if (p_valid) {
          if (p->cov >= conf->mincov) {
            fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%1.3f\n", p->chrm, p->pos-1, p->end, p->beta);
          }
          p_valid = 0;
        }

        if (!merged && parse_ok && (et&ET_HCG)) {
          bed1_t *tmp;
          tmp = p; p = b; b = tmp;
          p_valid=1;
        }
      }
      str.l = 0;                /* clean line */
      if (c==EOF) break;
    } else {
      kputc(c, &str);
    }
  }
  if (p_valid) {
    if (p->cov >= conf->mincov) {
      fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%1.3f\n", p->chrm, p->pos-1, p->end, p->beta);
    }
  }

}

void vcf2gch(gzFile FH, conf_t *conf) {

  bed1_t *b=init_bed1();
  kstring_t str;
  str.s = 0; str.l = str.m = 0;
  uint8_t et;
  while (1) {
    int c=gzgetc(FH);
    if (c=='\n' || c==EOF) {
      if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) {
        int parse_ok = vcf_parse1(str.s, b, &et);
        if (parse_ok && (et&ET_GCH) && b->cov >= conf->mincov) {
          fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%1.3f\n", b->chrm, b->pos-1, b->end, b->beta);
        }
      }
      str.l = 0;                /* clean line */
      if (c==EOF) break;
    } else {
      kputc(c, &str);
    }
  }
}

void vcf2c(gzFile FH, conf_t *conf) {

  bed1_t *b=init_bed1();
  kstring_t str;
  str.s = 0; str.l = str.m = 0;
  uint8_t et;
  while (1) {
    int c=gzgetc(FH);
    if (c=='\n' || c==EOF) {
      if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) {
        int parse_ok = vcf_parse1(str.s, b, &et); 
        if (parse_ok && (et&ET_C) && b->cov >= conf->mincov) {
          fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%1.3f\n", b->chrm, b->pos-1, b->end, b->beta);
        }
      }
      str.l = 0;                /* clean line */
      if (c==EOF) break;
    } else {
      kputc(c, &str);
    }
  }
}


int main_vcf2bed(int argc, char *argv[]) { 
  conf_t conf = {.verbose=0, .mincov=3};
  strcpy(conf.target, "cg");

  int c;
  while ((c = getopt(argc, argv, "k:t:V:h")) >= 0) {
    switch (c) {
    case 'k': conf.mincov = atoi(optarg); break;
    case 't': {
      if (strcmp(optarg, "cg")  !=0 && 
          strcmp(optarg, "hcg") !=0 && 
          strcmp(optarg, "gch") !=0 && 
          strcmp(optarg, "c")   !=0   ) {
        fprintf(stderr, "[%s:%d] Invalid option for -t: %s. \n", __func__, __LINE__, optarg);
        fflush(stderr);
        exit(1);
      }
      strcpy(conf.target, optarg);
      break;
    }
    case 'V': conf.verbose = atoi(optarg); break;
    case 'h': {
      fprintf(stderr, "\n");
      fprintf(stderr, "Usage: biscuit vcf2bed [options] vcf \n");
      fprintf(stderr, "Input options:\n");
      fprintf(stderr, "     -t STRING extract type {cg, hcg, gch, c} [%s]\n", conf.target);
      fprintf(stderr, "     -k INT    minimum coverage [%d]\n", conf.mincov);
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
    fprintf(stderr, "[%s:%d] Please provide input vcf\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  char *vcf_fn = argv[optind];
  gzFile FH=gzopen(vcf_fn, "r");
  if (!FH) {
    fprintf(stderr, "[%s:%d] Cannot open VCF file: %s\n", __func__, __LINE__, vcf_fn);
    fflush(stderr);
    exit(1);
  }
  if (strcmp(conf.target, "cg")==0) vcf2cg(FH, &conf);
  if (strcmp(conf.target, "hcg")==0) vcf2hcg(FH, &conf);
  if (strcmp(conf.target, "gch")==0) vcf2gch(FH, &conf);
  if (strcmp(conf.target, "c")==0) vcf2c(FH, &conf);

  return 0;
}
