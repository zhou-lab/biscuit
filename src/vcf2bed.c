#include <stdlib.h>
#include <string.h>
#include <zlib.h>
/* #include "wstr.h" */
#include "wzbed.h"
#include "wzvcf.h"


typedef struct conf_t {
  char target[5];
  int mincov;
  int showcov;
  int splitcg;
  int verbose;
  int destrand;
} conf_t;

typedef struct bed_vcfdata_t {
  char ref;
  int nsamples;
  double *betas;
  int *covs;
  char *cx;              /* context */
} bed_vcfdata_t;

static void free_bed_data(void *data) {
  bed_vcfdata_t *bd = (bed_vcfdata_t*) data;
  free(bd->cx);
  free(bd->betas);
  free(bd->covs);
  free(data);
}

static bed1_t *init_bed1(int nsamples) {
  bed1_t *b = init_bed1_core(NULL);
  b->data = calloc(1, sizeof(bed_vcfdata_t));
  bed_vcfdata_t *bd = (bed_vcfdata_t*) b->data;
  bd->nsamples = nsamples;
  bd->betas = calloc(bd->nsamples, sizeof(double));
  bd->covs = calloc(bd->nsamples, sizeof(int));
  return b;
}

#define free_bed1(b) free_bed1_core(b, free_bed_data)

/* typedef struct bed1_t { */
/*   char *chrm; */
/*   int64_t pos; */
/*   int64_t end; */
/*   char ref; */
/*   double beta; */
/*   int cov; */

/*   /\* reference and */
/*    * alternative as in VCF *\/ */
/*   char *vcfref; */
/*   char *vcfalt; */
/* } bed1_t; */

/* bed1_t *init_bed1() { */
/*   bed1_t *b=calloc(1, sizeof(bed1_t)); */
/*   b->chrm = 0; */
/*   b->vcfref = 0; */
/*   b->vcfalt = 0; */
/*   return b; */
/* } */

/* void free_bed1(bed1_t *b) { */
/*   free(b); */
/*   free(b->chrm); */
/*   free(b->vcfref); */
/*   free(b->vcfalt); */
/* } */

static int pass_coverage(bed1_t *b, conf_t *conf) {
  bed_vcfdata_t *bd = (bed_vcfdata_t*) b->data;
  int i;
  for (i=0; i<bd->nsamples; ++i)
    if (bd->covs[i] >= conf->mincov)
      return 1;
  return 0;
}

/* output one cytosine */
static void format_cytosine_bed1(bed1_t *b, conf_t *conf, target_v *targets) {

  if (b == NULL || b->tid < 0) return;
  if (!pass_coverage(b, conf)) return;
  
  fprintf(stdout, "%s\t%"PRId64"\t%"PRId64, target_name(targets, b->tid), b->beg, b->end);
  int i;
  bed_vcfdata_t *bd = (bed_vcfdata_t*) b->data;
  for (i=0; i<bd->nsamples; ++i) {
    if (bd->betas[i] < 0)
      fputs("\t.", stdout);
    else
      fprintf(stdout, "\t%1.3f", bd->betas[i]);
    if (conf->showcov) {
      fprintf(stdout, "\t%d", bd->covs[i]);
    }
  }
  fprintf(stdout, "\t%c\t%s", bd->ref, bd->cx);
  putchar('\n');
}

const char *genotypes[] = {"0/0", "0/1", "1/1"};

/* output SNP */
static void format_snp_bed1(bed1_t *b, conf_t *conf, char *vcfref, char *vcfalt, target_v *targets, char **fmt_gt, char **fmt_sp, int *allele_sp, int n_alleles) {

  if (b == NULL || b->tid < 0) return;

  /* compute highest non-ref AF and coverage */
  int i; int cov=0, haltcnt=0;
  for (i=0; i<n_alleles; ++i) {
    cov += allele_sp[i];
    if (i)
      if (allele_sp[i] > haltcnt)
        haltcnt = allele_sp[i];
  }
  if (cov < conf->mincov) return;

  bed_vcfdata_t *bd = (bed_vcfdata_t*) b->data;
  fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%s\t%s", target_name(targets, b->tid), b->beg, b->end, vcfref, vcfalt);
  for (i = 0; i<bd->nsamples; ++i)
    fprintf(stdout, "\t%s\t%s\t%d\t%1.2f", fmt_gt[i], fmt_sp[i], cov, (double) haltcnt / cov);
  putchar('\n');
}

/* output C,G in a CpG context together */
static void format_cytosine_bed2(bed1_t *p, bed1_t *b, conf_t *conf, target_v *targets) {

  if (p != NULL && p->tid < 0) p = NULL;
  if (b != NULL && b->tid < 0) b = NULL;
  if (b == NULL && p == NULL) return;

  bed_vcfdata_t *bd = b ? (bed_vcfdata_t*) b->data : NULL;
  bed_vcfdata_t *pd = p ? (bed_vcfdata_t*) p->data : NULL;

  bed1_t *tmp;
  if ((pd != NULL && pd->ref == 'G') || (bd != NULL && bd->ref == 'C')) {
    tmp = b; b = p; p = tmp;
    bd = b ? (bed_vcfdata_t*) b->data : NULL;
    pd = p ? (bed_vcfdata_t*) p->data : NULL;
  }

  if (pd != NULL && pd->ref != 'C') return;
  if (bd != NULL && bd->ref != 'G') return;
  if (pd != NULL && bd != NULL && pd->nsamples != bd->nsamples) {
    wzfatal("[%s:%d] Error: sample number does't match (%d vs %d).", __func__, __LINE__, pd->nsamples, bd->nsamples);
  }
  int nsamples = bd ? bd->nsamples : pd->nsamples;

  /* the merged coverage and betas */
  int *covs2 = calloc(nsamples, sizeof(int));
  double *betas2 = calloc(nsamples, sizeof(double));

  int i, cov; double ret;
  int cov_pass = 0;
  for (i=0; i<nsamples; ++i) {
    cov = 0; ret = 0;
    if (bd) {
      cov += bd->covs[i];
      ret += bd->covs[i] * bd->betas[i];
    }
    if (pd) {
      cov += pd->covs[i];
      ret += pd->covs[i] * pd->betas[i];
    }
    covs2[i] = cov;
    if (cov > 0) {
      betas2[i] = ret / (double) cov;
    } else {
      betas2[i] = -1;
    }
    if (cov >= conf->mincov) cov_pass = 1;
  }

  if (!cov_pass) goto CLEANUP;
  int tid = p ? p->tid : b->tid;
  int64_t beg = p ? p->beg : b->beg-1;
  int64_t end = b ? b->end : p->end+1;
  char *cx = pd ? pd->cx : bd->cx;
  fprintf(stdout, "%s\t%"PRId64"\t%"PRId64, target_name(targets, tid), beg, end);
  for (i=0; i<nsamples; ++i) {
    if (betas2[i] >= 0)
      fprintf(stdout, "\t%1.3f", betas2[i]);
    else
      fprintf(stdout, "\t.");
    if (conf->showcov) {
      fprintf(stdout, "\t%d", covs2[i]);

      /* display the C and G in separate columns */
      if (conf->splitcg) {
        if (pd && pd->covs[i] > 0) {
          fprintf(stdout, "\tC:%1.3f:%d", pd->betas[i], pd->covs[i]);
        } else {
          fprintf(stdout, "\tC:.:0");
        }
        if (bd && bd->covs[i] > 0) {
          fprintf(stdout, "\tG:%1.3f:%d", bd->betas[i], bd->covs[i]);
        } else {
          fprintf(stdout, "\tG:.:0");
        }
      }
    }
  }
  if (conf->showcov) {          /* show reference and sequence context */
    fprintf(stdout, "\tCG\t%s", cx);
  }
  putchar('\n');


 CLEANUP:
  free(covs2);
  free(betas2);
}


/* #define ET_C   0x1 */
/* #define ET_CG  0x2 */
/* #define ET_HCG 0x4 */
/* #define ET_GCH 0x8 */

/* static int vcf_parse1(char *line, bed1_t *b, uint8_t *et, char *cx) { */

/*   char *tok; */
/*   char *linerest=0, *fieldrest=0; */

/*   /\* CHROM *\/ */
/*   tok=strtok_r(line, "\t", &linerest); */
/*   b->chrm = realloc(b->chrm, strlen(tok)+1); */
/*   strcpy(b->chrm, tok); */

/*   /\* POS *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   ensure_number(tok); */
/*   b->pos = atoi(tok); b->end = b->pos; */
  
/*   /\* ID *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */

/*   /\* REF *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   b->ref = toupper(tok[0]); */
/*   b->vcfref = realloc(b->vcfref, strlen(tok)+1); */
/*   strcpy(b->vcfref, tok); */

/*   /\* ALT *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   b->vcfalt = realloc(b->vcfalt, strlen(tok)+1); */
/*   strcpy(b->vcfalt, tok); */

/*   /\* QUAL *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */

/*   /\* PASS FILTER *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */

/*   /\* INFO *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   *et=0; */
/*   int i; */
/*   for (i=0; i<(signed)strlen(tok)-2; ++i) { */
/*     if (strncmp(tok+i,"CX=",3)==0) { */
/*       memcpy(cx, tok+i+3, 3); */
/*       if (cx[2]==';' || cx[2]=='\t') cx[2]='\0'; */
/*       else cx[3]='\0'; */
/*     } */
/*     if (strncmp(tok+i,"N5=",3)==0) { */
/*       if (tok[i+5]=='C') { */
/*         *et |= ET_C; */
/*         if (tok[i+4]=='G' && tok[i+6]!='G') *et |= ET_GCH; */
/*         if (tok[i+6]=='G') { */
/*           *et |= ET_CG; */
/*           if (tok[i+4]!='G') *et |= ET_HCG; */
/*         } */
/*       } */
/*       break; */
/*     } */
/*   } */

/*   /\* not any of the target *\/ */
/*   if (!*et) return 0; */

/*   /\* FORMAT *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   int cv_index=-1, bt_index=-1; */
/*   char *field; */
/*   field=strtok_r(tok, ":", &fieldrest); */
/*   i=0; */
/*   while(field) { */
/*     if (field[0]=='C' && field[1]=='V') { */
/*       cv_index = i; */
/*     } */
/*     if (field[0]=='B' && field[1]=='T') { */
/*       bt_index = i; */
/*     } */
/*     ++i; */
/*     field=strtok_r(NULL, ":", &fieldrest); */
/*   } */

/*   /\* no coverage or beta info *\/ */
/*   if (cv_index<0 || bt_index<0) return 0; */

/*   /\* FORMAT content *\/ */
/*   tok = strtok_r(NULL, "\t", &linerest); */
/*   i=0; */
/*   int coverage=-1; double beta=-1.0; */
/*   field=strtok_r(tok, ":", &fieldrest); */
/*   while (field) { */
/*     if (i==cv_index) { */
/*       ensure_number(field); */
/*       coverage=atoi(field); */
/*     } */
/*     if (i==bt_index) { */
/*       ensure_number(field); */
/*       beta=atof(field); */
/*     } */
/*     ++i; */
/*     field=strtok_r(NULL, ":", &fieldrest); */
/*   } */

/*   /\* no coverage or beta info *\/ */
/*   if (coverage<0||beta<0) return 0; */

/*   b->cov = coverage; */
/*   b->beta = beta; */

/*   return 1; */
/* } */

/* Convert vcf_record_t to bed1_t */
static void vcf_record2bed1(bed1_t *b, vcf_record_t *rec, vcf_file_t *vcf) {

  char *info_cx = get_vcf_record_info("CX", rec->info);

  char **fmt_bt; int n_fmt_bt;
  char **fmt_cv; int n_fmt_cv;
  get_vcf_record_fmt("BT", rec->fmt, vcf, &fmt_bt, &n_fmt_bt);
  get_vcf_record_fmt("CV", rec->fmt, vcf, &fmt_cv, &n_fmt_cv);

  if (fmt_bt != NULL && n_fmt_bt != vcf->n_tsamples)
    wzfatal("[%s:%d] Invalid VCF file.", __func__, __LINE__);

  if (fmt_cv != NULL && n_fmt_cv != vcf->n_tsamples)
    wzfatal("[%s:%d] Invalid VCF file.", __func__, __LINE__);

  b->tid = rec->tid;
  b->beg = rec->pos-1;
  b->end = rec->pos;
  bed_vcfdata_t *bd = (bed_vcfdata_t*) b->data;
  bd->ref = rec->ref[0];
  bd->nsamples = vcf->n_tsamples;
  int i;

  if (n_fmt_bt == 0) {
    for (i=0; i<bd->nsamples; ++i)
      bd->betas[i] = -1.0;
  } else {
    for (i=0; i<n_fmt_bt; ++i) {
      if (!is_number(fmt_bt[i]) || strcmp(fmt_bt[i], ".") == 0)
        bd->betas[i] = -1.0;
      else
        bd->betas[i] = atof(fmt_bt[i]);
    }
  }

  if (n_fmt_cv == 0) {
    for (i=0; i<bd->nsamples; ++i)
      bd->covs[i] = 0;
  } else {
    for (i=0; i<n_fmt_cv; ++i) {
      if (!is_number(fmt_cv[i]) || strcmp(fmt_cv[i], ".") == 0)
        bd->covs[i] = 0;
      else
        bd->covs[i] = atoi(fmt_cv[i]);
    }
  }
  
  if (info_cx != NULL)
    bd->cx = strcpy_realloc(bd->cx, info_cx);
  else {
    free(bd->cx);
    bd->cx = NULL;
  }

  free(info_cx);
  free_char_array(fmt_bt, n_fmt_bt);
  free_char_array(fmt_cv, n_fmt_cv);
}

/* cx can be either "CG" or "HCG" */
static void vcf2cg(vcf_file_t *vcf, conf_t *conf, const char *cx) {

  bed1_t *b=init_bed1(vcf->n_tsamples);        /* base */
  bed1_t *p=init_bed1(vcf->n_tsamples);        /* previous base */
  b->tid = -1;
  p->tid = -1;

  vcf_record_t *rec = init_vcf_record();
  bed_vcfdata_t *bd, *pd;
  while (vcf_read_record(vcf, rec)) {

    vcf_record2bed1(b, rec, vcf);
    bd = (bed_vcfdata_t*) b->data;
    if (bd->cx == NULL || strcmp(bd->cx, cx) != 0) {
      b->tid = -1;
      continue;
    }
    
    if (conf->destrand) {
      bd = (bed_vcfdata_t*) b->data;
      pd = (bed_vcfdata_t*) p->data;
      if (p->tid >= 0 && b->tid >= 0 &&
          p->tid == b->tid && p->beg+1 == b->beg &&
          pd->ref == 'C' && bd->ref == 'G') { /* report C,G together */
        format_cytosine_bed2(p, b, conf, vcf->targets);
        b->tid = -1;  /* invalidate current, to avoid double report */
      } else { /* report previous base alone */
        format_cytosine_bed2(p, NULL, conf, vcf->targets);
      }
    } else {
      format_cytosine_bed1(p, conf, vcf->targets);
    }

    bed1_t *tmp;
    tmp = p; p = b; b = tmp;
  }

  if (conf->destrand)
    format_cytosine_bed2(p, NULL, conf, vcf->targets);
  else
    format_cytosine_bed1(p, conf, vcf->targets);

  free_bed1(b);
  free_bed1(p);
  free_vcf_record(rec);
}


void vcf2gch(vcf_file_t *vcf, conf_t *conf) {

  bed1_t *b=init_bed1(vcf->n_tsamples);
  vcf_record_t *rec = init_vcf_record();
  while (vcf_read_record(vcf, rec)) {
    vcf_record2bed1(b, rec, vcf);
    bed_vcfdata_t *bd = (bed_vcfdata_t*) b->data;
    if (bd->cx == NULL || strcmp(bd->cx, "GCH") != 0)
      continue;
    format_cytosine_bed1(b, conf, vcf->targets);
  }

  free_bed1(b);
  free_vcf_record(rec);
}

static void vcf2c(vcf_file_t *vcf, conf_t *conf) {

  bed1_t *b=init_bed1(vcf->n_tsamples);
  vcf_record_t *rec = init_vcf_record();
  while (vcf_read_record(vcf, rec)) {
    vcf_record2bed1(b, rec, vcf);
    bed_vcfdata_t *bd = (bed_vcfdata_t*) b->data;
    if (bd->ref != 'C' && bd->ref != 'G')
      continue;
    format_cytosine_bed1(b, conf, vcf->targets);
  }

  free_bed1(b);
  free_vcf_record(rec);
}


static void vcf2snp(vcf_file_t *vcf, conf_t *conf) {

  bed1_t *b=init_bed1(vcf->n_tsamples);
  vcf_record_t *rec = init_vcf_record();
  while (vcf_read_record(vcf, rec)) {
    vcf_record2bed1(b, rec, vcf);

    if (strcmp(rec->alt, ".") != 0) {
      char **fmt_gt; int n_fmt_gt;
      char **fmt_sp; int n_fmt_sp;
      get_vcf_record_fmt("GT", rec->fmt, vcf, &fmt_gt, &n_fmt_gt);
      get_vcf_record_fmt("SP", rec->fmt, vcf, &fmt_sp, &n_fmt_sp);

      /* parse out all alleles */
      char **alleles; int n_alleles;
      line_get_fields(rec->alt, ",", &alleles, &n_alleles);
      
      alleles = realloc(alleles, (++n_alleles)*sizeof(char*));
      memmove(alleles + 1, alleles, (n_alleles-1)*sizeof(char*));
      alleles[0] = strdup(rec->ref);

      /* parse out allele support in each sample */
      int *allele_sp = calloc(n_alleles*n_fmt_sp, sizeof(int));
      int j;
      for (j=0; j<n_fmt_sp; ++j) {
        int n_allele_sppairs; char **allele_sppairs;
        line_get_fields(fmt_sp[j], ",", &allele_sppairs, &n_allele_sppairs);
        int k;
        for (k=0; k<n_allele_sppairs; ++k) {
          char *ae;               /* pointing to first digit of allele count */
          for (ae = allele_sppairs[k]; !isdigit(*ae); ++ae);
          int ai;
          for (ai=0; ai<n_alleles; ++ai)
            if (strncmp(alleles[ai], allele_sppairs[k], ae-allele_sppairs[k]) == 0)
              break;
          if (ai < n_alleles)
            allele_sp[j*n_alleles+ai] = atoi(ae);
          else {
            wzfatal("Allele %s not found in %s\n", allele_sppairs[k], vcf->line);
          }
        }
        free_char_array(allele_sppairs, n_allele_sppairs);
      }

      format_snp_bed1(b, conf, rec->ref, rec->alt, vcf->targets, fmt_gt, fmt_sp, allele_sp, n_alleles);
      free(allele_sp);
      free_char_array(alleles, n_alleles);
      free_char_array(fmt_gt, n_fmt_gt);
      free_char_array(fmt_sp, n_fmt_sp);
    }
  }

  free_bed1(b);
  free_vcf_record(rec);
}

/* void vcf2ch(gzFile FH, char *target_samples, conf_t *conf) { */

/*   bed1_t *b=init_bed1(); */
/*   kstring_t str; */
/*   str.s = 0; str.l = str.m = 0; */
/*   uint8_t et; */
/*   char cx[4] = ""; */
/*   while (1) { */
/*     int c=gzgetc(FH); */
/*     if (c=='\n' || c==EOF) { */
/*       if (str.l>2 && str.s[0] == '#' str.s */
/*       if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) { */
/*         int parse_ok = vcf_parse1(str.s, b, &et, cx); */
/*         if (parse_ok && (et&ET_C) && !(et&ET_CG) && b->cov >= conf->mincov) */
/*           format_cytosine_bed1(b, conf, cx); */
/*       } */
/*       str.l = 0;                /\* clean line *\/ */
/*       if (c==EOF) break; */
/*     } else { */
/*       kputc(c, &str); */
/*     } */
/*   } */
/* } */

/* void vcf2gch(gzFile FH, conf_t *conf) { */

/*   bed1_t *b=init_bed1(); */
/*   kstring_t str; */
/*   str.s = 0; str.l = str.m = 0; */
/*   uint8_t et; */
/*   char cx[4] = ""; */
/*   while (1) { */
/*     int c=gzgetc(FH); */
/*     if (c=='\n' || c==EOF) { */
/*       if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) { */
/*         int parse_ok = vcf_parse1(str.s, b, &et, cx); */
/*         if (parse_ok && (et&ET_GCH) && b->cov >= conf->mincov) */
/*           format_cytosine_bed1(b, conf, cx); */
/*       } */
/*       str.l = 0;                /\* clean line *\/ */
/*       if (c==EOF) break; */
/*     } else { */
/*       kputc(c, &str); */
/*     } */
/*   } */
/* } */

/* void vcf2c(gzFile FH, conf_t *conf) { */

/*   bed1_t *b=init_bed1(); */
/*   kstring_t str; */
/*   str.s = 0; str.l = str.m = 0; */
/*   uint8_t et; */
/*   char cx[4] = ""; */
/*   while (1) { */
/*     int c=gzgetc(FH); */
/*     if (c=='\n' || c==EOF) { */
/*       if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) { */
/*         int parse_ok = vcf_parse1(str.s, b, &et, cx);  */
/*         if (parse_ok && (et&ET_C) && b->cov >= conf->mincov) */
/*           format_cytosine_bed1(b, conf, cx); */
/*       } */
/*       str.l = 0;                /\* clean line *\/ */
/*       if (c==EOF) break; */
/*     } else { */
/*       kputc(c, &str); */
/*     } */
/*   } */
/* } */

/* static int vcf_parse_snp1(char *line, bed1_t *b) { */

/*   int i; */
/*   char *tok; */
/*   char *linerest=0, *fieldrest=0; */
  
/*   /\* CHROM *\/ */
/*   tok=strtok_r(line, "\t", &linerest); */
/*   b->chrm = realloc(b->chrm, strlen(tok)+1); */
/*   strcpy(b->chrm, tok); */

/*   /\* POS *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   ensure_number(tok); */
/*   b->pos = atoi(tok); b->end = b->pos; */
  
/*   /\* ID *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */

/*   /\* REF *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   b->ref = toupper(tok[0]); */
/*   b->vcfref = realloc(b->vcfref, strlen(tok)+1); */
/*   strcpy(b->vcfref, tok); */

/*   if (strlen(b->vcfref) != 1) return 0; */
  
/*   /\* ALT *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   b->vcfalt = realloc(b->vcfalt, strlen(tok)+1); */
/*   strcpy(b->vcfalt, tok); */

/*   if (strlen(b->vcfalt) != 1) return 0; */

/*   /\* QUAL *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */

/*   /\* PASS FILTER *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */

/*   if (strcmp(tok, "PASS") != 0) return 0; */

/*   /\* INFO *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */

/*   /\* FORMAT *\/ */
/*   tok=strtok_r(NULL, "\t", &linerest); */
/*   int dp_index=-1, gt_index=-1; */
/*   char *field; */
/*   field=strtok_r(tok, ":", &fieldrest); */
/*   i=0; */
/*   while(field) { */
/*     if (field[0]=='D' && field[1]=='P') { */
/*       dp_index = i; */
/*     } */
/*     if (field[0]=='G' && field[1]=='T') { */
/*       gt_index = i; */
/*     } */
/*     ++i; */
/*     field=strtok_r(NULL, ":", &fieldrest); */
/*   } */

/*   /\* no coverage or beta info *\/ */
/*   if (dp_index<0 || gt_index<0) return 0; */

/*   /\* FORMAT content *\/ */
/*   tok = strtok_r(NULL, "\t", &linerest); */
/*   i=0; */
/*   int coverage=-1; int is_var=0; */
/*   field=strtok_r(tok, ":", &fieldrest); */
/*   while (field) { */
/*     if (i==dp_index) { */
/*       ensure_number(field); */
/*       coverage=atoi(field); */
/*     } */
/*     if (i==gt_index) { */
/*       if (strchr(field, '1')!=NULL) */
/*         is_var = 1; */
/*     } */
/*     ++i; */
/*     field=strtok_r(NULL, ":", &fieldrest); */
/*   } */

/*   /\* no coverage or beta info *\/ */
/*   if (coverage<0 || !is_var) return 0; */

/*   b->cov = coverage; */

/*   return 1; */
/* } */


/* void vcf2snp(gzFile FH, conf_t *conf) { */
  
/*   bed1_t *b=init_bed1(); */
/*   kstring_t str; */
/*   str.s = 0; str.l = str.m = 0; */
/*   while (1) { */
/*     int c=gzgetc(FH); */
/*     if (c=='\n' || c==EOF) { */
/*       if (str.l>2 && str.s[0] != '#' && strcount_char(str.s, '\t')>=8) { */
/*         int parse_ok = vcf_parse_snp1(str.s, b); */
/*         if (parse_ok && b->cov >= conf->mincov) */
/*           format_snp_bed1(b, conf); */
/*       } */
/*       str.l = 0;                /\* clean line *\/ */
/*       if (c==EOF) break; */
/*     } else { */
/*       kputc(c, &str); */
/*     } */
/*   } */
/*   free_bed1(b); */
/* } */

int main_vcf2bed(int argc, char *argv[]) { 
  conf_t conf = {.verbose=0, .mincov=3, .destrand=1, .showcov=0, .splitcg = 0};
  strcpy(conf.target, "cg");
  char *target_samples = NULL;

  int c;
  while ((c = getopt(argc, argv, "k:t:cupV:s:h")) >= 0) {
    switch (c) {
    case 'k': conf.mincov = atoi(optarg); break;
    case 't': {
      if (strcmp(optarg, "c") !=0 &&
          strcmp(optarg, "cg") !=0 &&
          strcmp(optarg, "hcg") !=0 && 
          strcmp(optarg, "gch") !=0 &&
          strcmp(optarg, "snp") !=0) {
        fprintf(stderr, "[%s:%d] Invalid option for -t: %s. \n", __func__, __LINE__, optarg);
        fflush(stderr);
        exit(1);
      }
      strcpy(conf.target, optarg);
      break;
    }
    case 'c': conf.showcov = 1; break;
    case 'u': conf.destrand = 0; break;
    case 'p': conf.splitcg = 1; break;
    case 'V': conf.verbose = atoi(optarg); break;
    case 's': target_samples = strdup(optarg); break;
    case 'h': {
      fprintf(stderr, "\n");
      fprintf(stderr, "Usage: biscuit vcf2bed [options] vcf \n");
      fprintf(stderr, "Input options:\n");
      fprintf(stderr, "     -t STRING extract type {c, cg, ch, hcg, gch, snp} [%s]\n", conf.target);
      fprintf(stderr, "     -k INT    minimum coverage [%d]\n", conf.mincov);
      fprintf(stderr, "     -s STRING sample, (takes \"FIRST\", \"LAST\", \"ALL\", or specific sample names separated by \",\")[FIRST]\n");
      fprintf(stderr, "     -c        show coverage, reference base(s) and sequence context as extra columns\n");
      fprintf(stderr, "     -u INT    suppress merging C and G in the CpG context (destrand, for cg and hcg, when both strands are in the hcg context).\n");
      fprintf(stderr, "     -p        display CpG together and show C and G as additional columns (reference_base:beta:depth).\n");
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

  if (!target_samples) {
    target_samples = strdup("FIRST");
  }

  if (optind >= argc) {
    fprintf(stderr, "[%s:%d] Please provide input vcf\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  char *vcf_fn = argv[optind];
  vcf_file_t *vcf = init_vcf_file(vcf_fn);
  index_vcf_samples(vcf, target_samples);

  if (strcmp(conf.target, "cg")==0) vcf2cg(vcf, &conf, "CG");
  if (strcmp(conf.target, "hcg")==0) vcf2cg(vcf, &conf, "HCG");
  if (strcmp(conf.target, "gch")==0) vcf2gch(vcf, &conf);
  if (strcmp(conf.target, "c")==0) vcf2c(vcf, &conf);
  if (strcmp(conf.target, "snp")==0) vcf2snp(vcf, &conf);

  free_vcf_file(vcf);
  free(target_samples);

  return 0;
}
