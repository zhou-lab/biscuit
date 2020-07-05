/* extract bed from VCF file
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
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
 * 
**/

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <zlib.h>
#include "wzbed.h"
#include "wzvcf.h"


typedef struct conf_t {
  char target[5];
  int mincov;
  int showctxt;
} conf_t;

typedef struct bed_data_t {
  char ref;
  int nsamples;
  double *betas;
  int *covs;
  char *cx;        // context
  char n5[6];      // 5 nucleotide context
} bed_data_t;

static void free_bed_data(void *data) {
  bed_data_t *bd = (bed_data_t*) data;
  free(bd->cx);
  free(bd->betas);
  free(bd->covs);
  free(data);
}

static void init_bed_data(bed1_t *b, void *aux_data) {
  b->data = calloc(1, sizeof(bed_data_t));
  bed_data_t *bd = (bed_data_t*) b->data;
  bd->nsamples = *((int*) aux_data);
  bd->betas = calloc(bd->nsamples, sizeof(double));
  bd->covs = calloc(bd->nsamples, sizeof(int));
}

static int pass_coverage(bed1_t *b, conf_t *conf) {
  bed_data_t *bd = (bed_data_t*) b->data;
  int i;
  for (i=0; i<bd->nsamples; ++i)
    if (bd->covs[i] >= conf->mincov)
      return 1;
  return 0;
}

const char *genotypes[] = {"0/0", "0/1", "1/1"};

/* Convert vcf_record_t to bed1_t */
static void vcf_record2bed1(bed1_t *b, vcf_record_t *rec, vcf_file_t *vcf) {

  char *info_cx = get_vcf_record_info("CX", rec->info);
  char *info_n5 = get_vcf_record_info("N5", rec->info);

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
  bed_data_t *bd = (bed_data_t*) b->data;
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

  if (info_n5 == NULL || strlen(info_n5) != 5) strcpy(bd->n5, "NNNNN");
  else strncpy(bd->n5, info_n5, 6);

  free(info_cx);
  free(info_n5);
  free_char_array(fmt_bt, n_fmt_bt);
  free_char_array(fmt_cv, n_fmt_cv);
}

static void vcf2bed_ctxt(vcf_file_t *vcf, conf_t *conf, const char *cx) {

  bed1_t *b=init_bed1(init_bed_data, &vcf->n_tsamples);
  vcf_record_t *rec = init_vcf_record();
  while (vcf_read_record(vcf, rec)) {
    vcf_record2bed1(b, rec, vcf);
    bed_data_t *bd = (bed_data_t*) b->data;

    // skip missing context
    if (bd->cx == NULL) continue;

    if (strcmp(cx, "C") == 0) { // targeting all C
      if (bd->ref != 'C' && bd->ref != 'G') continue;
    } else if (strcmp(cx, "CH") == 0) { // targeting CH
      if (strcmp(bd->cx, "CHH") != 0 && strcmp(bd->cx, "CHG") != 0) continue;
    } else if (strcmp(bd->cx, cx) != 0) continue; // all other cases

    if (b == NULL || b->tid < 0) continue;
    if (!pass_coverage(b, conf)) continue;

    // chrm, beg, end
    fprintf(stdout, "%s\t%"PRId64"\t%"PRId64, target_name(vcf->targets, b->tid), b->beg, b->end);
    // ctxt
    if (conf->showctxt) fprintf(stdout, "\t%c\t%s\t%.2s\t%.5s", bd->ref, bd->cx?bd->cx:"NA", bd->n5+2, bd->n5);
    int i;
    for (i=0; i<bd->nsamples; ++i) {
      // betas
      if (bd->betas[i] < 0) fputs("\t.", stdout);
      else fprintf(stdout, "\t%1.3f", bd->betas[i]);
      // coverage
      fprintf(stdout, "\t%d", bd->covs[i]);
    }
    if (fputc('\n', stdout) < 0 && errno == EPIPE) exit(1);
  }

  free_bed1(b, free_bed_data);
  free_vcf_record(rec);
}

static void vcf2bed_snp(vcf_file_t *vcf, conf_t *conf) {

  bed1_t *b=init_bed1(init_bed_data, &vcf->n_tsamples);
  vcf_record_t *rec = init_vcf_record();
  while (vcf_read_record(vcf, rec)) {
    vcf_record2bed1(b, rec, vcf);

    bed_data_t *bd = (bed_data_t*) b->data;
    if (strcmp(rec->alt, ".") != 0) {
      char **fmt_gt; int n_fmt_gt;
      char **fmt_sp; int n_fmt_sp;
      char **fmt_ac; int n_fmt_ac;
      char **fmt_af; int n_fmt_af;
      get_vcf_record_fmt("GT", rec->fmt, vcf, &fmt_gt, &n_fmt_gt);
      get_vcf_record_fmt("SP", rec->fmt, vcf, &fmt_sp, &n_fmt_sp);
      get_vcf_record_fmt("AC", rec->fmt, vcf, &fmt_ac, &n_fmt_ac);
      get_vcf_record_fmt("AF1", rec->fmt, vcf, &fmt_af, &n_fmt_af);

      if (n_fmt_sp != bd->nsamples || n_fmt_gt != bd->nsamples || n_fmt_ac != bd->nsamples || n_fmt_af != bd->nsamples)
        wzfatal("Malformed VCF file (unmatched no. records) in %s\n", vcf->line);

      /* parse out all alleles */
      // char **alleles; int n_alleles;
      // use the AB tag for alt
      // char *AB = get_vcf_record_info("AB", rec->info);
      /* line_get_fields(AB, ",", &alleles, &n_alleles); */
      /* free(AB); */
      /*        */
      /* alleles = realloc(alleles, (++n_alleles)*sizeof(char*)); */
      /* memmove(alleles + 1, alleles, (n_alleles-1)*sizeof(char*)); */
      /* alleles[0] = strdup(rec->ref); */
      /*  */
      /* [> parse out allele support in each sample <] */
      /* int j; */
      /* int *allele_sp = calloc(n_alleles*n_fmt_sp, sizeof(int)); */
      /* for (j=0; j<n_fmt_sp; ++j) { */
      /*   if (strcmp(fmt_sp[j], ".")==0) continue; */
      /*   int n_allele_sppairs; char **allele_sppairs; */
      /*   line_get_fields(fmt_sp[j], ",", &allele_sppairs, &n_allele_sppairs); */
      /*   int k; */
      /*   for (k=0; k<n_allele_sppairs; ++k) { */
      /*     // pointing to first digit of allele count */
      /*     char *ae; for (ae = allele_sppairs[k]; !isdigit(*ae); ++ae); */
      /*  */
      /*     // locate which allele */
      /*     int ai; */
      /*     for (ai=0; ai<n_alleles; ++ai) */
      /*       if (strncmp(alleles[ai], allele_sppairs[k], ae-allele_sppairs[k]) == 0) break; */
      /*      */
      /*     if (ai < n_alleles) */
      /*       allele_sp[j*n_alleles+ai] = atoi(ae); */
      /*     else */
      /*       wzfatal("Allele %s not found in %s\n", allele_sppairs[k], vcf->line); */
      /*   } */
      /*   free_char_array(allele_sppairs, n_allele_sppairs); */
      /* } */
      /*  */
      if (b == NULL || b->tid < 0) goto END;

      /* compute highest non-ref AF and coverage */
      int highest_cov = 0; int sid;
      double highest_af = 0.0;
      for (sid=0; sid<bd->nsamples; ++sid) {
        int cov = atoi(fmt_ac[sid]);
        if (cov > highest_cov) highest_cov = cov;
        double af = atof(fmt_af[sid]); // upon failure atof gives 0.0, which is OK
        if (af > highest_af) highest_af = af;
      }
      if (highest_cov < conf->mincov) goto END;
      if (highest_af <= 0.0) goto END;

      // output
      // chrm, beg, end, ref, alt
      fprintf(stdout, "%s\t%"PRId64"\t%"PRId64"\t%s\t%s", 
          target_name(vcf->targets, b->tid), b->beg, b->end, rec->ref, rec->alt);
      // genotype, support, cov, vaf
      for (sid = 0; sid<bd->nsamples; ++sid) {
        /* int highest_altcnt=0, cov=0; // recomputed cov, could be more efficient here */
        /* int *allele_sp1 = allele_sp + sid*n_alleles; */
        /* for (i=0; i<n_alleles; ++i) { */
        /*   if (i && allele_sp[i] > highest_altcnt) highest_altcnt = allele_sp1[i]; */
        /*   cov += allele_sp1[i]; */
        /* } */
        putchar('\t'); fputs(fmt_gt[sid], stdout);
        putchar('\t'); fputs(fmt_sp[sid], stdout);
        putchar('\t'); fputs(fmt_ac[sid], stdout);
        putchar('\t'); fputs(fmt_af[sid], stdout);
        /* if (cov) fprintf(stdout, "%1.2f", (double) highest_altcnt / cov); */
        /* else putchar('.'); */
      }

      if (fputc('\n', stdout) < 0 && errno == EPIPE) exit(1);

      // putchar('\n');

      // free(allele_sp);
      // free_char_array(alleles, n_alleles);
END:
      free_char_array(fmt_gt, n_fmt_gt);
      free_char_array(fmt_sp, n_fmt_sp);
      free_char_array(fmt_ac, n_fmt_ac);
      free_char_array(fmt_af, n_fmt_af);
    }
  }

  free_bed1(b, free_bed_data);
  free_vcf_record(rec);
}

static int usage(conf_t *conf) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit vcf2bed [options] <in.vcf>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -t STR    Extract type {c, cg, ch, hcg, gch, snp} [%s]\n", conf->target);
    fprintf(stderr, "    -k INT    Minimum coverage [%d]\n", conf->mincov);
    fprintf(stderr, "    -s STR    Sample, (takes \"FIRST\", \"LAST\", \"ALL\", or specific\n");
    fprintf(stderr, "                  sample names separated by \",\") [FIRST]\n");
    fprintf(stderr, "    -e        Show context (reference base, context group {CG,CHG,CHH},\n");
    fprintf(stderr, "                  2-base {CA,CC,CG,CT} and 5-base context) before beta\n");
    fprintf(stderr, "                  value and coverage column\n");
    fprintf(stderr, "    -h        This help\n");
    fprintf(stderr, "\n");

    return 1;
}

int main_vcf2bed(int argc, char *argv[]) { 
  conf_t conf = {.mincov=3, .showctxt=0};
  strcpy(conf.target, "CG");
  char *target_samples = NULL;

  int c;
  if (argc<2) return usage(&conf);
  while ((c = getopt(argc, argv, ":t:k:s:eh")) >= 0) {
      switch (c) {
          case 'k': conf.mincov = atoi(optarg); break;
          case 't': {
              if (strlen(optarg) > 4) wzfatal("Invalid option for -t: %s.\n", optarg);
              strcpy(conf.target, optarg);
              break;
          }
          case 's': target_samples = strdup(optarg); break;
          case 'e': conf.showctxt = 1; break;
          case 'h': return usage(&conf); break;
          case ':': usage(&conf); wzfatal("Option needs an argument: -%c\n", optopt);
          case '?': usage(&conf); wzfatal("Unrecognized option: -%c\n", optopt);
          default: usage(&conf); wzfatal("Unrecognized option: -%c\n",c); break;
      }
  }

  // default to only FIRST sample
  if (!target_samples) target_samples = strdup("FIRST");

  if (optind >= argc) { usage(&conf); wzfatal("Please provide input vcf.\n"); }
  vcf_file_t *vcf = init_vcf_file(argv[optind]);
  index_vcf_samples(vcf, target_samples);

  char *raw_target = strdup(conf.target);
  wzstrupr(conf.target); // use uppercase internally
  if (strcasecmp(conf.target, "CG") != 0 &&
      strcasecmp(conf.target, "CH") != 0 &&
      strcasecmp(conf.target, "C") != 0 &&
      strcasecmp(conf.target, "HCG") != 0 &&
      strcasecmp(conf.target, "GCH") != 0 &&
      strcasecmp(conf.target, "SNP") != 0) wzfatal("Invalid option for -t: %s.\n", raw_target);
  free(raw_target);

  if (strcmp(conf.target, "SNP")==0) vcf2bed_snp(vcf, &conf);
  else vcf2bed_ctxt(vcf, &conf, conf.target);

  free_vcf_file(vcf);
  free(target_samples);

  return 0;
}
