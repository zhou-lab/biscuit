#include "pileup.h"

/* typedef enum {BSS_RETENTION, BSS_CONVERSION, BSS_OTHER} bsstate_t; */
/* typedef enum {MCT, MCG, MCA, MGT, MGC, MGA} mutation_t; */
/* const char alts[] = "TGATCA"; */
const char nt256int8_to_mutcode[6] = "ACGTYR";
const char *cytosine_context[] = {"CG","CHG","CHH","CN"};

typedef struct {
  int n;                        /* number of sites */
  pileup_data_v **data;
} pileup_t;

pileup_t *init_pileup(int n) {
  pileup_t *p = malloc(sizeof(pileup_t));
  p->n = n;
  p->data = calloc(n, sizeof(pileup_data_v*));
  return p;
}

void destroy_pileup(pileup_t *p) {
  int i;
  for (i=0; i<p->n; ++i) {
    if (p->data[i]) free_pileup_data_v(p->data[i]);
  }
  free(p->data);
  free(p);
}

typedef struct {
  char *bam_fn;                 /* on stack */
  char *ref_fn;                 /* on stack */
  wqueue_t(window) *q;
  wqueue_t(record) *rq;
  conf_t *conf;
} result_t;

void pop_record_by_block_id(record_v *records, int64_t block_id, record_t *record) {
  uint64_t i;
  record_t *r;
  for (i=0; i<records->size; ++i) {
    r = ref_record_v(records, i);
    if (r->block_id == block_id) {
      *record = *r;             /* copy the record and set slot on shelf to OBSOLETE */
      r->block_id = RECORD_SLOT_OBSOLETE;
      return;
    }
  }
  record->block_id = RECORD_SLOT_OBSOLETE;
}

void put_into_record_v(record_v *records, record_t rec) {
  uint64_t i;
  record_t *r;

  /* fill blanks */
  for (i=0; i<records->size; ++i) {
    r = ref_record_v(records, i);
    if (r->block_id == RECORD_SLOT_OBSOLETE) {
      *r = rec;
      return;
    }
  }

  /* get a new slot */
  r = next_ref_record_v(records);
  *r = rec;
  return;
}

void merge_bsrate(bsrate_t *t, bsrate_t *s) {

  int i;
  if (t->m==0) bsrate_init(t, s->m);
  for (i=0; i<t->m; ++i) {
    t->ct_unconv[i] += s->ct_unconv[i];
    t->ct_conv[i] += s->ct_conv[i];
    t->ga_unconv[i] += s->ga_unconv[i];
    t->ga_conv[i] += s->ga_conv[i];
    t->ct_unconv_m[i] += s->ct_unconv_m[i];
    t->ct_conv_m[i] += s->ct_conv_m[i];
    t->ga_unconv_m[i] += s->ga_unconv_m[i];
    t->ga_conv_m[i] += s->ga_conv_m[i];
  }
}

void *write_func(void *data) {
  writer_conf_t *c = (writer_conf_t*) data;

  FILE *out;
  if (c->outfn) out=fopen(c->outfn, "w");
  else out=stdout;

  FILE *stats;
  if (c->statsfn) {
    stats=fopen(c->statsfn, "w");
  } else if (c->outfn) {
    c->statsfn = calloc(strlen(c->outfn)+7,1);
    strcpy(c->statsfn, c->outfn);
    strcat(c->statsfn, ".stats");
    stats=fopen(c->statsfn, "w");
  } else {
    stats=stderr;
  }
  
  if (c->header) fputs(c->header, out);
  int64_t next_block = 0;
  record_v *records = init_record_v(20);

  /* statistics */
  int64_t *l=(int64_t*)calloc(c->targets->size, sizeof(int64_t));
  int64_t *n=(int64_t*)calloc(c->targets->size, sizeof(int64_t));
  int64_t *n_uniq=(int64_t*)calloc(c->targets->size, sizeof(int64_t));

  double *betasum_context=(double*)calloc(c->targets->size*3, sizeof(double));
  int64_t *cnt_context=(int64_t*)calloc(c->targets->size*3, sizeof(int64_t));
  int i;
  bsrate_t b;
  bsrate_init(&b, c->conf->bsrate_max_pos);

  while (1) {
    record_t rec;
    wqueue_get(record, c->q, &rec);
    if(rec.block_id == RECORD_QUEUE_END) break;
    if (rec.block_id == next_block) {
      do {
        if (rec.s.s) {
          fputs(rec.s.s, out);

          /* statistics */
          l[rec.tid] += rec.l;
          n[rec.tid] += rec.n;
          n_uniq[rec.tid] += rec.n_uniq;

          /* methlevelaverage */
          for (i=0; i<3; ++i) {
            betasum_context[rec.tid*3+i] += rec.betasum_context[i];
            cnt_context[rec.tid*3+i] += rec.cnt_context[i];
          }

          /* merge bsrate */
          merge_bsrate(&b, &rec.b);
        }
        free(rec.s.s);
        bsrate_free(&rec.b);

        /* get next block from shelf if available else return OBSOLETE 
           and retrieve new block from queue  */
        next_block++;
        pop_record_by_block_id(records, next_block, &rec);
      } while (rec.block_id != RECORD_SLOT_OBSOLETE);
    } else {                    /* shelf the block if not next */
      put_into_record_v(records, rec);
    }
  }

  /* output statistics */
  fprintf(stats, "\n#### base coverage ####\n");
  fprintf(stats, "chrom\tlen\tcov\tcov_uniq\n");
  uint32_t k; int64_t l0=0, n0=0, n_uniq0=0;
  for (k=0; k<c->targets->size; ++k) {
    if (l[k] == 0) continue;
    fprintf(stats, "%s\t%"PRId64"\t%"PRId64"\t%1.2f%%\t%"PRId64"\t%1.2f%%\n", get_target_v(c->targets, k).name, l[k], n[k], (double)n[k]/(double)l[k]*100, n_uniq[k], (double)n_uniq[k]/(double)l[k]*100);
    l0 += l[k]; n0 += n[k]; n_uniq0 += n_uniq[k];
  }
  fprintf(stats, "whole_genome\t%"PRId64"\t%"PRId64"\t%1.2f%%\t%"PRId64"\t%1.2f%%\n", l0, n0, (double)n0/(double)l0*100, n_uniq0, (double)n_uniq0/(double)l0*100);

  fprintf(stats, "\n#### methlevelaverage ####\n");
  fprintf(stats, "chrom\tCGn\tCGb\tCHGn\tCHGb\tCHHn\tCHHb\n");
  double betasum_context0[3] = {0.0,0.0,0.0};
  int64_t cnt_context0[3] = {0,0,0}; int k0;
  for (k=0; k<c->targets->size; ++k) {
    if (l[k] == 0) continue;    /* skip chrom with no base coverage */
    fprintf(stats, "%s\t%"PRId64"\t%1.3f%%\t%"PRId64"\t%1.3f%%\t%"PRId64"\t%1.3f%%\n",
            get_target_v(c->targets, k).name,
            cnt_context[k*3], betasum_context[k*3] / (double) cnt_context[k*3]*100,
            cnt_context[k*3+1], betasum_context[k*3+1] / (double) cnt_context[k*3+1]*100,
            cnt_context[k*3+2], betasum_context[k*3+2] / (double) cnt_context[k*3+2]*100);
    for (k0=0; k0<3; ++k0) {
      cnt_context0[k0] += cnt_context[k*3+k0];
      betasum_context0[k0] += betasum_context[k*3+k0];
    }
  }
  fprintf(stats, "whole_genome\t%"PRId64"\t%1.3f%%\t%"PRId64"\t%1.3f%%\t%"PRId64"\t%1.3f%%\n",
          cnt_context0[0], betasum_context0[0] / (double) cnt_context0[0]*100,
          cnt_context0[1], betasum_context0[1] / (double) cnt_context0[1]*100,
          cnt_context0[2], betasum_context0[2] / (double) cnt_context0[2]*100);

  fprintf(stats, "\n#### bisculfite conversion rate ####\n");
  fprintf(stats, "pos\tct_c\tct_u\tct_r\tga_c\tga_u\tga_r\tct_cm\tct_um\tct_rm\tga_cm\tga_um\tga_rm\n");
  for (i=0; i<b.m; ++i){
    if (!b.ct_conv[i] && !b.ga_conv[i] && !b.ct_conv_m[i] && !b.ga_conv_m[i]) continue;
    fprintf(stats, "%d\t%d\t%d\t%1.3f%%\t%d\t%d\t%1.3f%%\t%d\t%d\t%1.3f%%\t%d\t%d\t%1.3f%%\n",
            i+1, b.ct_conv[i], b.ct_unconv[i],
            (double) b.ct_conv[i] / (double)(b.ct_conv[i]+b.ct_unconv[i])*100,
            b.ga_conv[i], b.ga_unconv[i],
            (double) b.ga_conv[i] / (double)(b.ga_conv[i]+b.ga_unconv[i])*100,
            b.ct_conv_m[i], b.ct_unconv_m[i],
            (double) b.ct_conv_m[i] / (double)(b.ct_conv_m[i]+b.ct_unconv_m[i])*100,
            b.ga_conv_m[i], b.ga_unconv_m[i],
            (double) b.ga_conv_m[i] / (double)(b.ga_conv_m[i]+b.ga_unconv_m[i])*100);
  }

  free(c->statsfn);
  free(l); free(n); free(n_uniq);
  free(cnt_context);
  free(betasum_context);
  bsrate_free(&b);
  
  free_record_v(records);
  if (c->outfn) {    /* for stdout, will close at the end of main */
    fflush(out);
    fclose(out);
  }
  return 0;
}


uint8_t infer_bsstrand(refseq_t *rs, bam1_t *b, uint32_t min_base_qual) {

  /* infer bsstrand from nC2T and nG2A on high quality bases */
  
  bam1_core_t *c = &b->core;
  uint32_t rpos = c->pos+1, qpos = 0;
  uint32_t op, oplen;
  char rb, qb;
  int i, nC2T=0, nG2A=0; unsigned j;
  for (i=0; i<c->n_cigar; ++i) {
    op = bam_cigar_op(bam1_cigar(b)[i]);
    oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
    switch(op) {
    case BAM_CMATCH:
      for (j=0; j<oplen; ++j) {
        rb = toupper(getbase_refseq(rs, rpos+j));
        qb = bscall(b, qpos+j);
        if (bam1_qual(b)[qpos+j] < min_base_qual) continue;
        if (rb == 'C' && qb == 'T') nC2T++;
        if (rb == 'G' && qb == 'A') nG2A++;
      }
      rpos += oplen;
      qpos += oplen;
      break;
    case BAM_CINS:
      qpos += oplen;
      break;
    case BAM_CDEL:
      rpos += oplen;
      break;
    case BAM_CSOFT_CLIP:
      qpos += oplen;
      break;
    case BAM_CHARD_CLIP:
      qpos += oplen;
      break;
    default:
      fprintf(stderr, "Unknown cigar, %u\n", op);
      abort();
    }
  }
  if (nC2T >= nG2A) return 0;
  else return 1;
}

uint8_t get_bsstrand(refseq_t *rs, bam1_t *b, uint32_t min_base_qual) {
  uint8_t *s = bam_aux_get(b, "ZS");
  if (s) {
    s++;
    if (*s == '+') return 0;
    else if (*s == '-') return 1;
  }

  s = bam_aux_get(b, "YD");     /* bwa-meth flag */
  if (s) {
    s++;
    if (*s == 'f') return 0;
    else if (*s == 'r') return 1;
  }

  /* otherwise, guess the bsstrand from nCT and nGA */
  return infer_bsstrand(rs, b, min_base_qual);
}

static void verbose_format(uint8_t bsstrand, pileup_data_v *dv, kstring_t *s) {

  uint32_t i, nf;

  /* return if no record match the bsstrand */
  int n=0;
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) ++n;
  }
  if (!n) return;

  char b='0'+bsstrand;

  /* 1. base */
  ksprintf(s, ";Bs%c=", b);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc(d->qb, s);
  }

  /* 2. status array */
  ksprintf(s, ";Sta%c=", b);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc('0'+d->stat, s);
  }

  /* 3. base quality */
  ksprintf(s, ";Bq%c=", b);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc(d->qual+33, s);
  }

  /* 4. strand */
  ksprintf(s, ";Str%c=", b);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc(d->strand?'-':'+', s);
  }

  /* 5. position on read */
  ksprintf(s, ";Pos%c=", b);
  nf = 0;
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) {
      if (nf) kputc(',', s);
      else nf = 1;
      kputuw(d->qpos, s);
    }
  }

  /* 6. retention count 
     retention count, for diagnosing incomplete converted
     reads from CpH sites and mitochondrial sites */
  ksprintf(s, ";Rret%c=", b);
  nf = 0;
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) {
      if (nf) kputc(',', s);
      else nf = 1;
      kputuw(d->cnt_ret, s);
    }
  }
}

cytosine_context_t fivenuc_context(refseq_t *rs, uint32_t rpos, char rb, char *fivenuc) {
  /* char fivenuc[5]; */
  if (rpos == 1) {
    subseq_refseq2(rs, 1, fivenuc+2, 3);
    fivenuc[0] = fivenuc[1] = 'N';
  } else if (rpos == 2) {
    subseq_refseq2(rs, 1, fivenuc+1, 4);
    fivenuc[0] = 'N';
  } else if (rpos == (unsigned) rs->seqlen) {
    subseq_refseq2(rs, rpos-2, fivenuc, 3);
    fivenuc[3] = fivenuc[4] = 'N';
  } else if (rpos == (unsigned) rs->seqlen-1) {
    subseq_refseq2(rs, rpos-2, fivenuc, 4);
    fivenuc[4] = 'N';
  } else {
    subseq_refseq2(rs, rpos-2, fivenuc, 5);
  }
  if (rb == 'G') nt256char_rev_ip(fivenuc, 5);

  if (fivenuc[3] == 'N') return CONTEXT_NA;
  else if (fivenuc[3] == 'G') return CONTEXT_CG;
  else if (fivenuc[4] == 'G') return CONTEXT_CHG;
  else if (fivenuc[4] == 'N') return CONTEXT_NA;
  else return CONTEXT_CHH;
}

void plp_getcnts(pileup_data_v *dv, conf_t *conf, int cnts[9], int *_cm1, int *_cm2) {

  if (!dv) {
    *_cm1 = -1; *_cm2 = -1;
    return;
  }

  uint32_t i;
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv, i);
    /* read-position-based filtering */
    if (d->qual < conf->min_base_qual) continue;
    if (d->qpos < conf->min_dist_end ||
        d->rlen < d->qpos + conf->min_dist_end) continue;
    cnts[d->stat]++;
  }

  /* reset ambiguous mutation if they can be disambiguated */
  if (cnts[BSS_MC]>0 || cnts[BSS_MT]>0) cnts[BSS_MY] = 0;
  if (cnts[BSS_MG]>0 || cnts[BSS_MA]>0) cnts[BSS_MR] = 0;

  /* pick the top 2 mutations */
  int cm1=-1, cm2=-1;
  for (i=0; i<6; i++) {
    if (cnts[i] > 0) {
      if (cm1<0) {
        cm1 = i;
      } else if (cnts[i]>cnts[cm1]) {
        cm2 = cm1;
        cm1 = i;
      } else if (cm2<0 || cnts[i]>cnts[cm2]) {
        cm2 = i;
      }
    }
  }
  *_cm1 = cm1; *_cm2 = cm2;
}

int reference_supp(int cnts[9]) {
  int cref = 0;
  if (cnts[BSS_N]) cref += cnts[BSS_N];
  if (cnts[BSS_RETENTION] || cnts[BSS_CONVERSION])
    cref += cnts[BSS_RETENTION] + cnts[BSS_CONVERSION];
  return cref;
}

void allele_supp(char rb, int cref, int cm1, int cm2, int cnts[9], kstring_t *s) {

  if (cref)
    ksprintf(s, "%c%d", rb, cref);

  if (cm1 >= 0) {
    if (cref) kputc(',',s);
    ksprintf(s, "%c%d", nt256int8_to_mutcode[cm1], cnts[cm1]);
    if (cm2 >= 0) {
      ksprintf(s, ",%c%d", nt256int8_to_mutcode[cm2], cnts[cm2]);
      int i;
      for (i=0; i<6; ++i) {
        if (cnts[i]>0 && i!= cm1 && i!= cm2) {
          ksprintf(s, ",%c%d", nt256int8_to_mutcode[i], cnts[i]);
        }
      }
    }
  }
}

void pileup_genotype(int cref, int altsupp, conf_t *conf, char gt[4], double *_gl0, double *_gl1, double *_gl2, double *_gq) {

  double gl0=-1, gl1=-1, gl2=-1, gq=-1;
  if (cref >=0 || altsupp >= 0) {
    gl0 = log(conf->prior0) + genotype_lnlik(HOMOREF, cref, altsupp, conf->error, conf->contam);
    gl1 = log(conf->prior1) + genotype_lnlik(HET, cref, altsupp, conf->error, conf->contam);
    gl2 = log(conf->prior2) + genotype_lnlik(HOMOVAR, cref, altsupp, conf->error, conf->contam);
    if (gl0>gl1) {
      if (gl0>gl2) {
        gq = pval2qual(1 - exp(gl0 - ln_sum3(gl0, gl1, gl2)));
        strcpy(gt, "0/0");
      } else {
        gq = pval2qual(1 - exp(gl2 - ln_sum3(gl0, gl1, gl2)));
        strcpy(gt, "1/1");
      }
    } else if (gl1>gl2) {
      gq = pval2qual(1 - exp(gl1 - ln_sum3(gl0, gl1, gl2)));
      strcpy(gt, "0/1");
    } else {
      gq = pval2qual(1 - exp(gl2 - ln_sum3(gl0, gl1, gl2)));
      strcpy(gt, "1/1");
    }
  }
  *_gl0 = gl0; *_gl1 = gl1; *_gl2 = gl2; *_gq = gq;
}

static void plp_format(refseq_t *rs, char *chrm, uint32_t rpos,
                       pileup_data_v *dv, conf_t *conf, record_t *rec) {

  kstring_t *s = &rec->s;
  uint32_t i;
  char rb = toupper(getbase_refseq(rs, rpos));

  int cnts[9] = {0};
  int cm1, cm2;
  plp_getcnts(dv, conf, cnts, &cm1, &cm2);

  /* if not SNP but no signal for BSS_RETENTION or BSS_CONVERSION,
     skip the print when in non-verbose mode */
  if (cm1 < 0 && !conf->verbose
      && cnts[BSS_RETENTION] == 0
      && cnts[BSS_CONVERSION] == 0) return;

  /* MY and MR do not interfere */
  uint8_t methcallable=0;
  if (cnts[BSS_RETENTION] + cnts[BSS_CONVERSION] > 0) {
    if (cnts[BSS_MT]==0 && rb == 'C') methcallable = 1;
    if (cnts[BSS_MA]==0 && rb == 'G') methcallable = 1;
  }

  /* CHROM, POS, ID, REF */
  ksprintf(s, "%s\t%u\t.\t%c\t", chrm, rpos, rb);

  /* ALT
     if BSW shows G->A or BSC shows C->T, then a SNP, 
     no methylation information is inferred */
  if (cm1 >= 0) {
    uint32_t supp[6];
    for (i=0; i<6; ++i) supp[i] = (cnts[i]<<4) | i;
    qsort(supp, 6, sizeof(uint32_t), compare_supp);
    int fst = 1;
    for (i=0; i<6; ++i) {
      if ((supp[i]>>4) > 0) {
        char m = nt256int8_to_mutcode[supp[i]&0xf];
        if (m != rb) {
          if (!fst) kputc(',', s);
          kputc(m, s);
          fst = 0;
        }
      }
    }
  } else {
    kputc('.', s);
  }

  char gt[4] = "./.";
  int cref = reference_supp(cnts);
  int altsupp = cm1 >= 0?cnts[cm1]:0;
  double gl0, gl1, gl2, gq=-1;
  pileup_genotype(cref, altsupp, conf, gt, &gl0, &gl1, &gl2, &gq);

  /* QUAL */
  ksprintf(s, "\t%1.2f", gq);
  if (gq > 1) {
    kputs("\tPASS\t", s);
  } else {
    kputs("\tLowQual\t", s);
  }

  /* INFO tags */
  cytosine_context_t ctt=CONTEXT_NA;
  kputs("NS=1", s);
  if (methcallable) {
    char fivenuc[5];
    ctt = fivenuc_context(rs, rpos, rb, fivenuc);
    ksprintf(s, ";CX=%s", cytosine_context[ctt]);
    ksprintf(s, ";N5=%.5s", fivenuc);
  }

  /* FORMAT */
  kputs("\tDP:GT:GP:GQ", s);
  if (cref || cm1 >=0) kputs(":SP", s);
  if (methcallable) kputs(":CV:BT", s);

  /* BY SAMPLE FORMAT CONTEXT */
  if (dv) ksprintf(s, "\t%u", dv->size);
  else kputs("\t0", s);

  if (gq>0) {
    ksprintf(s, ":%s:%1.0f,%1.0f,%1.0f:%1.0f", gt, min(1000, -gl0), min(1000, -gl1), min(1000, -gl2), gq);
  } else {
    ksprintf(s, ":./.:.:.");
  }

  if (cref || cm1>=0) {
    kputc(':', s);
    allele_supp(rb, cref, cm1, cm2, cnts, s);
  }
  
  if (methcallable) {
    double beta = (double) cnts[BSS_RETENTION] / (double) (cnts[BSS_RETENTION]+cnts[BSS_CONVERSION]);
    if (ctt != CONTEXT_NA) {
      rec->betasum_context[ctt] += beta;
      rec->cnt_context[ctt]++;
    }
    ksprintf(s, ":%d:%1.2f", cnts[BSS_RETENTION]+cnts[BSS_CONVERSION], beta);
  }

  /* additional information printed on verbose theoretically this should be 
     put to FORMAT since they are sample-specific, but in FORMAT is usually 
     harder to parse, so they are appended to INFO these are not intended for
     formal submission, just for diagnostic purposes.
  */
  if (conf->verbose) {
    kputs("\tDIAGNOSE", s);
    if (methcallable)
      ksprintf(s, ";RN=%d;CN=%d", cnts[BSS_RETENTION], cnts[BSS_CONVERSION]);
    verbose_format(0, dv, s);
    verbose_format(1, dv, s);
  }

  kputc('\n', s);
}

/* return -1 if abnormal (missing bsstrand) */
uint32_t cnt_retention(refseq_t *rs, bam1_t *b, uint8_t bsstrand) {
  uint32_t cnt = 0;

  bam1_core_t *c = &b->core;
  uint32_t rpos = c->pos+1, qpos = 0;
  uint32_t op, oplen;
  char rb, qb;
  int i; unsigned j;
  for (i=0; i<c->n_cigar; ++i) {
    op = bam_cigar_op(bam1_cigar(b)[i]);
    oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
    switch(op) {
    case BAM_CMATCH:
      for (j=0; j<oplen; ++j) {
        rb = toupper(getbase_refseq(rs, rpos+j));
        qb = bscall(b, qpos+j);
        if (bsstrand) {
          if (rb == 'C' && qb == 'C') cnt++;
        } else {
          if (rb == 'G' && qb == 'G') cnt++;
        }
      }
      rpos += oplen;
      qpos += oplen;
      break;
    case BAM_CINS:
      qpos += oplen;
      break;
    case BAM_CDEL:
      rpos += oplen;
      break;
    case BAM_CSOFT_CLIP:
      qpos += oplen;
      break;
    case BAM_CHARD_CLIP:
      qpos += oplen;
      break;
    default:
      fprintf(stderr, "Unknown cigar, %u\n", op);
      abort();
    }
  }

  return cnt;
}

static void read_update_basecov(bam1_t *b, int *basecov, int *basecov_uniq, uint32_t wbeg, uint32_t wend, conf_t *conf) {

  bam1_core_t *c = &b->core;
  int i; unsigned j;
  uint32_t rpos = c->pos+1, qpos = 0;
  for (i=0; i<c->n_cigar; ++i) {
    uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
    uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
    switch(op) {
    case BAM_CMATCH:
      for (j=0; j<oplen; ++j) {
        if (rpos+j<wbeg || rpos+j>=wend) continue; /* include begin but not end */
        if (c->qual >= conf->min_mapq) basecov_uniq[rpos+j-wbeg]++;
        basecov[rpos+j-wbeg]++;
      }
      rpos += oplen;
      qpos += oplen;
      break;
    case BAM_CINS:
      qpos += oplen;
      break;
    case BAM_CDEL:
      rpos += oplen;
      break;
    case BAM_CSOFT_CLIP:
      qpos += oplen;
      break;
    case BAM_CHARD_CLIP:
      qpos += oplen;
      break;
    default:
      fprintf(stderr, "Unknown cigar, %u\n", op);
      abort();
    }
  }
}

static void calc_bsrate(uint8_t bsstrand, char rb, char qb, refseq_t *rs, uint32_t qp, uint32_t rp, uint8_t is_mito, bam1_core_t *c, bsrate_t *b, conf_t *conf) {

  int pos_on_template;
  if (c->flag & BAM_FREAD2) {
    if (abs(c->isize) > (signed) qp) pos_on_template = abs(c->isize) - qp;
    else return;
  } else pos_on_template = qp;

  if (pos_on_template >= conf->bsrate_max_pos) return;

  if (bsstrand && rb == 'G') {
    if (qb == 'A') {
      if (is_mito) { /* mitochondrial */
        b->ga_conv_m[pos_on_template]++;
      } else if (rp-1 > rs->beg && toupper(getbase_refseq(rs, rp-1)) != 'C') { /* non-mitochondrial, look at CpH context */
        b->ga_conv[pos_on_template]++;
      }
    } else if (qp == 'G') {
      if (is_mito) {          /* mitochondrial */
        b->ga_unconv_m[pos_on_template]++;
      } else if (rp-1 > rs->beg && toupper(getbase_refseq(rs, rp-1)) != 'C') { /* non-mitochondrial, look at CpH context */
        b->ga_unconv[pos_on_template]++;
      }
    }
  } else if (!bsstrand && rb == 'C') {
    if (qb == 'T') {
      if (is_mito) {          /* mitochondrial */
        b->ct_conv_m[pos_on_template]++;
      } else if (rp+1 < rs->end && toupper(getbase_refseq(rs, rp+1)) != 'G') { /* non-mitochondrial, look at CpH context */
        b->ct_conv[pos_on_template]++;
      }
    } else if (qp == 'C') {
      if (is_mito) {          /* mitochondrial */
        b->ct_unconv_m[pos_on_template]++;
      } else if (rp+1 < rs->end && toupper(getbase_refseq(rs, rp+1)) != 'G') { /* non-mitochondrial, look at CpH context */
        b->ct_unconv[pos_on_template]++;
      }
    }
  }
}

static void *process_func(void *data) {

  result_t *res = (result_t*) data;
  conf_t *conf = (conf_t*) res->conf;
  samfile_t *in = samopen(res->bam_fn, "rb", 0);
  bam_index_t *idx = bam_index_load(res->bam_fn);
  refseq_t *rs = init_refseq(res->ref_fn, 1000, 1000);
  int i; uint32_t j;

  record_t rec;
  window_t w;
  while (1) {

    wqueue_get(window, res->q, &w);
    if (w.tid == -1) break;

    rec.tid = w.tid;
    bsrate_init(&rec.b, conf->bsrate_max_pos);
    for (i=0; i<3; ++i) {
      rec.betasum_context[i] = 0;
      rec.cnt_context[i] = 0;
    }

    pileup_t *plp = init_pileup(w.end - w.beg);

    int *basecov = calloc(w.end-w.beg, sizeof(int));
    int *basecov_uniq = calloc(w.end-w.beg, sizeof(int));

    char *chrm = in->header->target_name[w.tid];
    uint8_t is_mito=0;
    if (strcmp(chrm, "chrM")==0 || strcmp(chrm, "MT")==0) is_mito=1;

    fetch_refseq(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
    bam_iter_t iter = bam_iter_query(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
    bam1_t *b = bam_init1();
    int ret;
    char qb, rb;
    while ((ret = bam_iter_read(in->x.bam, iter, b))>0) {

      /* uint8_t *bsstrand = bam_aux_get(b, "ZS"); */
      /* if (!bsstrand) continue; */
      /* bsstrand++; */

      uint8_t bsstrand = get_bsstrand(rs, b, conf->min_base_qual);
      read_update_basecov(b, basecov, basecov_uniq, w.beg, w.end, conf);
      
      /* read-based filtering */
      bam1_core_t *c = &b->core;
      if (c->qual < conf->min_mapq) continue;
      if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->min_read_len) continue;
      if (c->flag > 0){         /* only when any flag is set */
        if (conf->filter_secondary && c->flag & BAM_FSECONDARY) continue;
        if (conf->filter_duplicate && c->flag & BAM_FDUP) continue;
        if (conf->filter_ppair && !(c->flag & BAM_FPROPER_PAIR)) continue;
        if (conf->filter_qcfail && c->flag & BAM_FQCFAIL) continue;
      }

      uint8_t *nm = bam_aux_get(b, "NM");
      if (nm && bam_aux2i(nm)>conf->max_nm) continue;
      uint32_t cnt_ret = cnt_retention(rs, b, bsstrand);
      if (cnt_ret > conf->max_retention) continue;

      uint32_t rpos = c->pos+1, qpos = 0;
      for (i=0; i<c->n_cigar; ++i) {
        uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
        uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
        switch(op) {
        case BAM_CMATCH:
          for (j=0; j<oplen; ++j) {
            if (rpos+j<w.beg || rpos+j>=w.end) continue; /* include begin but not end */
            rb = toupper(getbase_refseq(rs, rpos+j));
            /* if (rb != 'C' && rb != 'G') continue; */
            qb = bscall(b, qpos+j);
            pileup_data_v **plp_data_vec = plp->data+rpos+j-w.beg;
            if (!*plp_data_vec) *plp_data_vec = init_pileup_data_v(2);
            pileup_data_t *d = next_ref_pileup_data_v(*plp_data_vec);
            d->qual = bam1_qual(b)[qpos+j];
            d->cnt_ret = (unsigned) cnt_ret;
            d->strand = (c->flag&BAM_FREVERSE)?1:0;
            d->qpos = qpos+j+1;
            d->rlen = c->l_qseq;
            d->bsstrand = bsstrand;
            d->qb = qb;

            calc_bsrate(bsstrand, rb, qb, rs, qpos+j, rpos+j, is_mito, c, &rec.b, conf);

            if (bsstrand) {          /* BSC */
              if (rb == 'G') {
                if (qb == 'A') d->stat = BSS_CONVERSION;
                else if (qb == 'G') d->stat = BSS_RETENTION;
                else d->stat = mutcode(qb);
              } else if (rb != qb) {
                if (qb == 'A') d->stat = BSS_MR;
                else d->stat = mutcode(qb);
              } else {
                d->stat = BSS_N;
              }
      
            } else {                    /* BSW */
              if (rb == 'C') {
                if (qb == 'T') d->stat = BSS_CONVERSION;
                else if (qb == 'C') d->stat = BSS_RETENTION;
                else d->stat = mutcode(qb);
              } else if (rb != qb) {
                if (qb == 'T') d->stat = BSS_MY;
                else d->stat = mutcode(qb);
              } else {
                d->stat = BSS_N;
              }
            }
          }
          rpos += oplen;
          qpos += oplen;
          break;
        case BAM_CINS:
          qpos += oplen;
          break;
        case BAM_CDEL:
          rpos += oplen;
          break;
        case BAM_CSOFT_CLIP:
          qpos += oplen;
          break;
        case BAM_CHARD_CLIP:
          qpos += oplen;
          break;
        default:
          fprintf(stderr, "Unknown cigar, %u\n", op);
          abort();
        }
      }
    }

    /* run through cytosines */
    rec.s.l = rec.s.m = 0; rec.s.s = 0; /* the record string */
    rec.block_id = w.block_id;
    for (j=w.beg; j<w.end; ++j) {
      rb = getbase_refseq(rs, j);
      pileup_data_v *plp_data = plp->data[j-w.beg];
      if (plp_data) {
        plp_format(rs, chrm, j, plp_data, res->conf, &rec);
      }
    }

    rec.l = w.end-w.beg;
    rec.n_uniq = 0; rec.n = 0;
    for (i=0; i<rec.l; ++i) {
      if (basecov[i] >= res->conf->min_cov) rec.n++;
      if (basecov_uniq[i] >= res->conf->min_cov) rec.n_uniq++;
    }

    /* put output string to output queue */
    wqueue_put2(record, res->rq, rec);

    destroy_pileup(plp);
    bam_destroy1(b);
    bam_iter_destroy(iter);
  }
  free_refseq(rs);
  samclose(in);
  bam_index_destroy(idx);
  return 0;
}

static void head_append_verbose(char *pb, char b, kstring_t *s) {
  ksprintf(s, "##FORMAT=<ID=Bs%c,Number=1,Type=String,Description=\"base identity, %s\">\n", b, pb);
  ksprintf(s, "##FORMAT=<ID=Sta%c,Number=1,Type=String,Description=\"Status code, %s (0,1,2,3 for mutation into A,C,G,T; 4,5 for Y,R; 6,7 for retention and conversion; 8 for other normal;)\">\n", b, pb);
  ksprintf(s, "##FORMAT=<ID=Bq%c,Number=1,Type=String,Description=\"base quality, %s\">\n", b, pb);
  ksprintf(s, "##FORMAT=<ID=Str%c,Number=1,Type=String;Description=\"strands, %s\">\n", b, pb);
  ksprintf(s, "##FORMAT=<ID=Pos%c,Number=1,Type=String;Description=\"position in read, %s\">\n", b, pb);
  ksprintf(s, "##FORMAT=<ID=Rret%c,Number=1,Type=String;Description=\"Number of retention in read, %s\">\n", b, pb);
}

static int usage(conf_t *conf) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: pileup [options] -r [ref.fa] -i [in.bam] -o [out.pileup] -g [chr1:123-234]\n");
  fprintf(stderr, "output format: chrm, pos, pos, refbase, mutation, cytosine_context, coverage, filtered_retention, filtered_conversion\n");
  fprintf(stderr, "Input options:\n\n");
  fprintf(stderr, "     -i        input bam.\n");
  fprintf(stderr, "     -r        reference in fasta.\n");
  fprintf(stderr, "     -g        region (optional, if not specified the whole bam will be processed).\n");
  fprintf(stderr, "     -s        step of window dispatching [%d].\n", conf->step);
  fprintf(stderr, "     -q        number of threads [%d] recommend 20.\n", conf->n_threads);
  fprintf(stderr, "\nOutputing format:\n\n");
  fprintf(stderr, "     -o        pileup output file [stdout]\n");
  fprintf(stderr, "     -w        pileup statistics output, e.g., bsrate, methlevelaverage etc. [stderr if no '-o' else [output].stats]\n");
  fprintf(stderr, "     -v        verbose (print additional info for diagnosis).\n");
  fprintf(stderr, "\nGenotyping parameters:\n\n");
  fprintf(stderr, "     -E        error rate [%1.3f].\n", conf->error);
  fprintf(stderr, "     -M        mutation rate [%1.3f].\n", conf->mu);
  fprintf(stderr, "     -C        contamination rate [%1.3f].\n", conf->contam);
  fprintf(stderr, "     -P        prior probability for heterozygous variant [%1.3f].\n", conf->prior1);
  fprintf(stderr, "     -Q        prior probability for homozygous variant [%1.3f].\n", conf->prior2);
  fprintf(stderr, "\nPileup filtering:\n\n");
  fprintf(stderr, "     -k        min read coverage in computing methlevelaverage and base coverage statistics [%d]\n", conf->min_cov);
  fprintf(stderr, "     -b        min base quality [%u].\n", conf->min_base_qual);
  fprintf(stderr, "     -m        minimum mapping quality [%u].\n", conf->min_mapq);
  fprintf(stderr, "     -t        max cytosine retention in a read [%u].\n", conf->max_retention);
  fprintf(stderr, "     -l        minimum read length [%u].\n", conf->min_read_len);
  fprintf(stderr, "     -e        minimum distance to end of a read [%u].\n", conf->min_dist_end);
  fprintf(stderr, "     -c        NO filtering secondary mapping.\n");
  fprintf(stderr, "     -u        NO filtering of duplicate.\n");
  fprintf(stderr, "     -p        NO filtering of improper pair (!BAM_FPROPER_PAIR).\n");
  fprintf(stderr, "     -S        bsrate maximum position [%d]\n", conf->bsrate_max_pos);
  fprintf(stderr, "     -n        maximum NM tag [%u].\n", conf->max_nm);
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");

  return 1;
}

void conf_init(conf_t *conf) {
  conf->step = 100000;
  conf->n_threads = 3;
  conf->bsrate_max_pos = 1000;
  conf->min_base_qual = 20;
  conf->min_mapq = 40;
  conf->min_cov = 3;
  conf->max_retention = 999999;
  conf->min_read_len = 10;
  conf->filter_qcfail = 1;
  conf->filter_secondary = 1;
  conf->filter_duplicate = 1;
  conf->filter_ppair = 1;
  conf->min_dist_end = 3;
  conf->max_nm = 255;
  conf->contam = 0.1;
  conf->error = 0.001;
  conf->mu = 0.001;
  conf->prior1 = 0.33333;
  conf->prior2 = 0.33333;
  conf->prior0 = 1.0 - conf->prior1 - conf->prior2;
  if (conf->prior0 < 0) {
    fprintf(stderr, "[Error] genotype prior0 (%1.3f) must be from 0 to 1. \n", conf->prior0);
    exit(1);
  }
  if (conf->prior1 < 0) {
    fprintf(stderr, "[Error] genotype prior1 (%1.3f) must be from 0 to 1. \n", conf->prior1);
    exit(1);
  }
  if (conf->prior2 < 0) {
    fprintf(stderr, "[Error] genotype prior2 (%1.3f) must be from 0 to 1. \n", conf->prior2);
    exit(1);
  }
  
  conf->verbose = 0;
  conf->noheader = 0;
}

int main_pileup(int argc, char *argv[]) {

  int c;
  char *reffn = 0;
  char *reg = 0;
  char *infn = 0;
  char *outfn = 0;
  char *statsfn = 0;
  conf_t conf;
  conf_init(&conf);

  if (argc<2) return usage(&conf);
  while ((c=getopt(argc, argv, "i:o:w:r:g:q:e:s:b:S:k:E:M:C:P:Q:t:n:m:l:cupvh"))>=0) {
    switch (c) {
    case 'i': infn = optarg; break;
    case 'o': outfn = optarg; break;
    case 'w': statsfn = strdup(optarg); break;
    case 'r': reffn = optarg; break;
    case 'g': reg = optarg; break;
    case 'q': conf.n_threads = atoi(optarg); break;
    case 's': conf.step = atoi(optarg); break;
    case 'e': conf.min_dist_end = atoi(optarg); break;
    case 'b': conf.min_base_qual = atoi(optarg); break;
    case 'S': conf.bsrate_max_pos = atoi(optarg); break;
    case 'k': conf.min_cov = atoi(optarg); break;
    case 'E': conf.error = atof(optarg); break;
    case 'M': conf.mu = atof(optarg); break;
    case 'C': conf.contam = atof(optarg); break;
    case 'P': conf.prior1 = atof(optarg); break;
    case 'Q': conf.prior2 = atof(optarg); break;
    case 't': conf.max_retention = atoi(optarg); break;
    case 'l': conf.min_read_len = atoi(optarg); break;
    case 'n': conf.max_nm = atoi(optarg); break;
    case 'm': conf.min_mapq = atoi(optarg); break;
    case 'c': conf.filter_secondary = 0; break;
    case 'u': conf.filter_duplicate = 0; break;
    case 'p': conf.filter_ppair = 0; break;
    case 'v': conf.verbose = 1; break;
    case 'h': return usage(&conf);
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
  }

  if (!infn || !reffn) {
    usage(&conf);
    exit(1);
  }

  wqueue_t(window) *wq = wqueue_init(window, 100000);
  pthread_t *processors = calloc(conf.n_threads, sizeof(pthread_t));
  result_t *results = calloc(conf.n_threads, sizeof(result_t));
  int i; unsigned j;
  samfile_t *in = samopen(infn, "rb", 0);

  /* process header */
  kstring_t header; header.l = header.m = 0; header.s = 0;
  kputs("##fileformat=VCFv4.1\n", &header);
  ksprintf(&header, "##reference=%s\n", reffn);
  ksprintf(&header, "##source=biscuitV%s\n", PACKAGE_VERSION);

  /* sort sequence name by alphabetic order, chr1, chr10, chr11 ... */
  target_v *targets = init_target_v(50);
  target_t *t;
  for (i=0; i<in->header->n_targets; ++i) {
    t = next_ref_target_v(targets);
    t->tid = i;
    t->name = in->header->target_name[i];
    t->len = in->header->target_len[i];
  }
  qsort(targets->buffer, targets->size, sizeof(target_t), compare_targets);
  for (j=0; j<targets->size; ++j) {
    t = ref_target_v(targets, j);
    ksprintf(&header, "##contig=<ID=%s,length=%d>\n", t->name, t->len);
  }
  kputs("##program=<cmd=biscuit", &header);
  for (i=0; i<argc; ++i)
    ksprintf(&header, " %s", argv[i]);
  kputs(">\n", &header);
  kputs("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n", &header);
  kputs("##INFO=<ID=CX,Number=1,Type=Float,Description=\"Cytosine context (CG, CHH or CHG)\">\n", &header);
  kputs("##INFO=<ID=N5,Number=1,Type=String,Description=\"5-nucleotide context, centered around target cytosine\">\n", &header);

  kputs("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n", &header);
  kputs("##FORMAT=<ID=SP,Number=R,Type=String,Description=\"Allele support (with filtering)\">\n", &header);
  kputs("##FORMAT=<ID=CV,Number=1,Type=Integer,Description=\"Effective (strand-specific) coverage on cytosine\">\n", &header);
  kputs("##FORMAT=<ID=BT,Number=1,Type=Float,Description=\"Cytosine methylation fraction (aka beta value, with filtering)\">\n", &header);
  kputs("##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype from normal\">\n", &header);
  kputs("##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods\">\n", &header);
  kputs("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality (phred-scaled)\">\n", &header);

  if (conf.verbose) {
    kputs("##FORMAT=<ID=RN,Number=1,Type=Integer,Description=\"Retention count (with filtering)\">\n", &header);
    kputs("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Conversion count (with filtering)\">\n", &header);
    char plpbsstrand[4];
    strcpy(plpbsstrand, "BSW");
    head_append_verbose(plpbsstrand, '0', &header);
    strcpy(plpbsstrand, "BSC");
    head_append_verbose(plpbsstrand, '1', &header);
  }

  kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n", &header);

  /* setup writer */
  pthread_t writer;
  writer_conf_t writer_conf = {
    .q = wqueue_init(record, 100000),
    .outfn = outfn,
    .statsfn = statsfn,
    .header = conf.noheader?0:header.s,
    .targets = targets,
    .conf = &conf,
  };
  pthread_create(&writer, NULL, write_func, &writer_conf);
  for (i=0; i<conf.n_threads; ++i) {
    results[i].q = wq;
    results[i].rq = writer_conf.q;
    results[i].ref_fn = reffn;
    results[i].bam_fn = infn;
    results[i].conf = &conf;
    pthread_create(&processors[i], NULL, process_func, &results[i]);
  }

  window_t w; memset(&w, 0, sizeof(window_t));
  uint32_t wbeg;
  int64_t block_id=0;

  /* process bam */
  if (reg) {                    /* regional */
    int tid;
    uint32_t beg, end;
    bam_parse_region(in->header, reg, &tid, (int*) &beg, (int*) &end);
    /* chromosome are assumed to be less than 2**29 */
    beg++; end++;
    if (beg<=0) beg = 1;
    if (end>in->header->target_len[tid]) end = in->header->target_len[tid];
    for (wbeg = beg; wbeg < end; wbeg += conf.step, block_id++) {
      w.tid = tid;
      w.block_id = block_id;
      w.beg = wbeg;
      w.end = wbeg + conf.step;
      if (w.end > end) w.end = end;
      wqueue_put(window, wq, &w);
    }
  } else {                      /* entire bam */
    for (j=0; j<targets->size; ++j) {
      t = ref_target_v(targets, j);
      for (wbeg = 1; wbeg < t->len; wbeg += conf.step, block_id++) {
        w.tid = t->tid;
        w.block_id = block_id;
        w.beg = wbeg;
        w.end = wbeg+conf.step;
        if (w.end > t->len) w.end = t->len;
        wqueue_put(window, wq, &w);
      }
    }
  }
  for (i=0; i<conf.n_threads; ++i) {
    w.tid = -1;
    wqueue_put(window, wq, &w);
  }

  for (i=0; i<conf.n_threads; ++i) {
    pthread_join(processors[i], NULL);
  }

  record_t rec = { .block_id = RECORD_QUEUE_END };
  wqueue_put2(record, writer_conf.q, rec);
  pthread_join(writer, NULL);
  wqueue_destroy(record, writer_conf.q);

  free_target_v(targets);
  free(results);
  free(processors);
  free(header.s);
  wqueue_destroy(window, wq);
  samclose(in);

  return 0;
}
