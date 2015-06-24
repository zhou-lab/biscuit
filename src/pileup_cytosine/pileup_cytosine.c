#include <ctype.h>
#include <stdlib.h>
#include "wqueue.h"
#include "encode.h"
#include "sam.h"
#include "refseq.h"
#include "kstring.h"
#include "wvec.h"

#define bscall(b, pos) bam_nt16_rev_table[bam1_seqi(bam1_seq(b), pos)]

typedef struct {
  int64_t block_id;
  int32_t tid;
  uint32_t beg, end;
} window_t;

DEFINE_WQUEUE(window, window_t)

typedef struct {
  int step;
  int n_threads;
  uint32_t min_base_qual;
  uint32_t max_retention;
  uint32_t min_read_len;
  uint8_t min_dist_end;
  uint8_t min_mapq;
  uint8_t max_nm;
  uint8_t filter_ppair:1;       /* filter BAM_FPROPER_PAIR */
  uint8_t filter_secondary:1;
  uint8_t filter_duplicate:1;
  uint8_t filter_qcfail:1;
  uint8_t verbose;
} conf_t;

/* typedef enum {BSS_RETENTION, BSS_CONVERSION, BSS_OTHER} bsstate_t; */
/* typedef enum {MCT, MCG, MCA, MGT, MGC, MGA} mutation_t; */
/* const char alts[] = "TGATCA"; */
const char nt256int8_to_mutcode[6] = "ACGTYR";
typedef enum {BSS_MA, BSS_MC, BSS_MG, BSS_MT,
              BSS_MY, BSS_MR, BSS_RETENTION, BSS_CONVERSION, BSS_N} status_t;

typedef struct {
  uint8_t bsstrand:1;
  uint8_t qual:7;
  uint8_t strand:1;
  uint16_t qpos;
  uint8_t cnt_ret;
  uint16_t rlen;                /* read length */
  char qb;
  status_t stat;                  /* code from mut-met status table */
} __attribute__((__packed__)) pileup_data_t;

DEFINE_VECTOR(pileup_data_v, pileup_data_t)

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

#define RECORD_QUEUE_END -2
#define RECORD_SLOT_OBSOLETE -1

typedef struct {
  int64_t block_id;
  kstring_t s;
} record_t;

DEFINE_VECTOR(record_v, record_t)

DEFINE_WQUEUE(record, record_t)

typedef struct {
  char *bam_fn;                 /* on stack */
  char *ref_fn;                 /* on stack */
  wqueue_t(window) *q;
  wqueue_t(record) *rq;
  conf_t *conf;
} result_t;

typedef struct {
  wqueue_t(record) *q;
  char *outfn;
} writer_conf_t;

void remove_record_by_block_id(record_v *records, int64_t block_id, record_t *record) {
  uint64_t i;
  record_t *r;
  for (i=0; i<records->size; ++i) {
    r = ref_record_v(records, i);
    if (r->block_id == block_id) {
      *record = *r;
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

void *write_func(void *data) {
  writer_conf_t *c = (writer_conf_t*) data;
  FILE *out;
  if (c->outfn) out=fopen(c->outfn, "w");
  else out=stdout;
  fflush(out);
  int64_t next_block = 0;
  record_v *records = init_record_v(20);
  while (1) {
    record_t rec;
    wqueue_get(record, c->q, &rec);
    if(rec.block_id == RECORD_QUEUE_END) break;
    if (rec.block_id == next_block) {
      do {
        if (rec.s.s) fputs(rec.s.s, out);
        free(rec.s.s);
        next_block++;
        remove_record_by_block_id(records, next_block, &rec);
      } while (rec.block_id != RECORD_SLOT_OBSOLETE);
    } else {
      put_into_record_v(records, rec);
    }
  }

  free_record_v(records);
  fclose(out);
  return 0;
}


uint8_t infer_bsstrand(refseq_t *rs, bam1_t *b, conf_t *conf) {

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
        if (bam1_qual(b)[qpos+j] < conf->min_base_qual) continue;
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

uint8_t get_bsstrand(refseq_t *rs, bam1_t *b, conf_t *conf) {
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
  return infer_bsstrand(rs, b, conf);
}

void verbose_format(uint8_t bsstrand, pileup_data_v *dv, kstring_t *s) {


  uint32_t i, nf;
  
  /* 1. base */
  kputc('\t', s);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc(d->qb, s);
  }

  /* 2. status array */
  kputc('\t', s);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc('0'+d->stat, s);
  }

  /* 3. base quality */
  kputc('\t', s);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc(d->qual+33, s);
  }

  /* 4. strand */
  kputc('\t', s);
  for (i=0; i<dv->size; ++i) {
    pileup_data_t *d = ref_pileup_data_v(dv,i);
    if (d->bsstrand == bsstrand) kputc(d->strand?'-':'+', s);
  }

  /* 5. position on read */
  kputc('\t', s);
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
  kputc('\t', s);
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

#define mutcode(a) (nt256char_to_nt256int8_table[(uint8_t)a])

void plp_format(refseq_t *rs, char *chrm, uint32_t rpos, pileup_data_v *dv, conf_t *conf, kstring_t *s) {
  uint32_t i;
  char rb = toupper(getbase_refseq(rs, rpos));

  int cnts[9] = {0};

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

  /* if not SNP but no signal for BSS_RETENTION or BSS_CONVERSION,
     skip the print when in non-verbose mode */
  if (cm1 < 0 && !conf->verbose
      && cnts[BSS_RETENTION] == 0
      && cnts[BSS_CONVERSION] == 0) return;

  /* MY and MR do not interfere */
  uint8_t methcallable=0;
  if (cnts[BSS_MT]==0 && rb == 'C') methcallable = 1;
  if (cnts[BSS_MA]==0 && rb == 'G') methcallable = 1;

  ksprintf(s, "%s\t%u\t%u\t%c", chrm, rpos-1, rpos, rb);
  
  /* if BSW shows G->A or BSC shows C->T, then a SNP, 
     no methylation information is inferred */
  if (cm1 >= 0) {
    kputc('\t', s);
    kputc(rb, s);
    kputc('>',s);
    int c = 0;
    if (cnts[BSS_N]) c += cnts[BSS_N];
    if (methcallable && (cnts[BSS_RETENTION] || cnts[BSS_CONVERSION]))
      c += cnts[BSS_RETENTION] + cnts[BSS_CONVERSION];
    if (c) ksprintf(s, "%c:%d,", rb, c);
    ksprintf(s, "%c:%d", nt256int8_to_mutcode[cm1], cnts[cm1]);
    if (cm2 >= 0) {
      ksprintf(s, ",%c:%d", nt256int8_to_mutcode[cm2], cnts[cm2]);
      for (i=0; i<6; ++i) {
        if (cnts[i]>0 && i!=(unsigned) cm1 && i!=(unsigned) cm2) {
          ksprintf(s, ",%c:%d", nt256int8_to_mutcode[i], cnts[i]);
        }
      }
    }
  } else {
    kputs("\t.", s);
  }
    
  if (methcallable) {
    /* for methylation output context */
    char trinuc[3];
    if (rpos == 1) {
      subseq_refseq2(rs, 1, trinuc+1, 2);
      trinuc[0] = 'N';
    } else if (rpos == (unsigned) rs->seqlen) {
      subseq_refseq2(rs, rpos-1, trinuc, 2);
      trinuc[2] = 'N';
    } else {
      subseq_refseq2(rs, rpos-1, trinuc, 3);
    }
    if (rb == 'G') {
      char trinuc_r[3];
      _nt256char_rev(trinuc_r, trinuc, 3);
      ksprintf(s, "\t%.3s", trinuc_r);
    } else {                    /* C,A,T context */
      ksprintf(s, "\t%.3s", trinuc);
    }
  } else {
    kputs("\t.", s);
  }
  
  /* coverage */
  ksprintf(s, "\t%u", dv->size);

  /* count retention and conversion */
  if (methcallable) {
    ksprintf(s, "\t%d\t%d", cnts[BSS_RETENTION], cnts[BSS_CONVERSION]);
  } else {
    kputs("\t.\t.", s);
  }

  /* additional information printed on verbose */
  if (conf->verbose) {
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


void *process_func(void *data) {

  result_t *res = (result_t*) data;
  conf_t *conf = (conf_t*) res->conf;
  samfile_t *in = samopen(res->bam_fn, "rb", 0);
  bam_index_t *idx = bam_index_load(res->bam_fn);
  refseq_t *rs = init_refseq(res->ref_fn, 1000, 1000);

  record_t rec;
  window_t w;
  int i; uint32_t j;
  while (1) {

    wqueue_get(window, res->q, &w);
    if (w.tid == -1) break;

    pileup_t *plp = init_pileup(w.end - w.beg);
    
    char *chrm = in->header->target_name[w.tid];
    fetch_refseq(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
    bam_iter_t iter = bam_iter_query(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
    bam1_t *b = bam_init1();
    int ret;
    char qb, rb;
    while ((ret = bam_iter_read(in->x.bam, iter, b))>0) {

      /* uint8_t *bsstrand = bam_aux_get(b, "ZS"); */
      /* if (!bsstrand) continue; */
      /* bsstrand++; */

      uint8_t bsstrand = get_bsstrand(rs, b, conf);

      bam1_core_t *c = &b->core;

      /* read-based filtering */
      if (c->qual < conf->min_mapq) continue;
      if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->min_read_len) continue;
      if (c->flag > 0){         /* only when any flag is set */
        if (conf->filter_secondary && c->flag & BAM_FSECONDARY) continue;
        if (conf->filter_duplicate && c->flag & BAM_FDUP) continue;
        if (conf->filter_ppair && !(c->flag & BAM_FPROPER_PAIR)) continue;
        if (conf->filter_qcfail && c->flag & BAM_FQCFAIL) continue;
      }
      
      uint32_t rpos = c->pos+1, qpos = 0;
      uint8_t *nm = bam_aux_get(b, "NM");
      if (nm && bam_aux2i(nm)>conf->max_nm) continue;
      uint32_t cnt_ret = cnt_retention(rs, b, bsstrand);
      if (cnt_ret > conf->max_retention) continue;

      rpos = c->pos+1; qpos = 0;
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
    rec.s.l = rec.s.m = 0; rec.s.s = 0;
    rec.block_id = w.block_id;
    for (j=w.beg; j<w.end; ++j) {
      rb = getbase_refseq(rs, j);
      pileup_data_v *plp_data = plp->data[j-w.beg];
      if (plp_data) {
        plp_format(rs, chrm, j, plp_data, res->conf, &rec.s);
      }
    }
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

typedef struct {
  int32_t tid;
  char *name;
  uint32_t len;
} target_t;

DEFINE_VECTOR(target_v, target_t)

int compare_targets(const void *a, const void *b) {
  return strcmp(((target_t*)a)->name, ((target_t*)b)->name);
}

int pileup_cytosine_main(char *reffn, char *infn, char *outfn, char *reg, conf_t *conf) {

  wqueue_t(window) *wq = wqueue_init(window, 100000);
  pthread_t *processors = calloc(conf->n_threads, sizeof(pthread_t));
  result_t *results = calloc(conf->n_threads, sizeof(result_t));
  int i; unsigned j;
  samfile_t *in = samopen(infn, "rb", 0);

  pthread_t writer;
  writer_conf_t writer_conf = {
    .q = wqueue_init(record, 100000),
    .outfn = outfn,
  };
  pthread_create(&writer, NULL, write_func, &writer_conf);
  for (i=0; i<conf->n_threads; ++i) {
    results[i].q = wq;
    results[i].rq = writer_conf.q;
    results[i].ref_fn = reffn;
    results[i].bam_fn = infn;
    results[i].conf = conf;
    pthread_create(&processors[i], NULL, process_func, &results[i]);
  }

  window_t w; memset(&w, 0, sizeof(window_t));
  uint32_t wbeg;
  int64_t block_id=0;
  target_v *targets = init_target_v(50);
  if (reg) {
    int tid;
    uint32_t beg, end;
    bam_parse_region(in->header, reg, &tid, (int*) &beg, (int*) &end);
    /* chromosome are assumed to be less than 2**29 */
    beg++; end++;
    if (beg<=0) beg = 1;
    if (end>in->header->target_len[tid]) end = in->header->target_len[tid];
    for (wbeg = beg; wbeg < end; wbeg += conf->step, block_id++) {
      w.tid = tid;
      w.block_id = block_id;
      w.beg = wbeg;
      w.end = wbeg + conf->step;
      if (w.end > end) w.end = end;
      wqueue_put(window, wq, &w);
    }
  } else {

    /* sort sequence name by alphabetic order, chr1, chr10, chr11 ... */
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
      for (wbeg = 1; wbeg < t->len; wbeg += conf->step, block_id++) {
        w.tid = t->tid;
        w.block_id = block_id;
        w.beg = wbeg;
        w.end = wbeg+conf->step;
        if (w.end > t->len) w.end = t->len;
        wqueue_put(window, wq, &w);
      }
    }
  }
  for (i=0; i<conf->n_threads; ++i) {
    w.tid = -1;
    wqueue_put(window, wq, &w);
  }

  for (i=0; i<conf->n_threads; ++i) {
    pthread_join(processors[i], NULL);
  }

  record_t rec = { .block_id = RECORD_QUEUE_END };
  wqueue_put2(record, writer_conf.q, rec);
  pthread_join(writer, NULL);
  wqueue_destroy(record, writer_conf.q);

  free_target_v(targets);
  free(results);
  free(processors);
  wqueue_destroy(window, wq);
  samclose(in);

  return 0;
}

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: pileup_cytosine [options] -r [ref.fa] -i [in.bam] -o [out.pileup] -g [chr1:123-234]\n");
  fprintf(stderr, "output format: chrm, pos, pos, refbase, mutation, cytosine_context, coverage, filtered_retention, filtered_conversion\n");
  fprintf(stderr, "Input options:\n");
  fprintf(stderr, "     -i        input bam.\n");
  fprintf(stderr, "     -r        reference in fasta.\n");
  fprintf(stderr, "     -g        region (optional, if not specified the whole bam will be processed).\n");
  fprintf(stderr, "     -o        pileup output file\n");
  fprintf(stderr, "     -s        step of window dispatching [100000].\n");
  fprintf(stderr, "     -q        number of threads [3] recommend 20.\n");
  fprintf(stderr, "     -b        min base quality [20].\n");
  fprintf(stderr, "     -m        minimum mapping quality [40].\n");
  fprintf(stderr, "     -t        max retention in a read [999999].\n");
  fprintf(stderr, "     -l        minimum read length [10].\n");
  fprintf(stderr, "     -e        minimum distance to end of a read [3].\n");
  fprintf(stderr, "     -c        NO filtering secondary mapping.\n");
  fprintf(stderr, "     -u        NO filtering of duplicate.\n");
  fprintf(stderr, "     -p        NO filtering of improper pair (!BAM_FPROPER_PAIR).\n");
  fprintf(stderr, "     -n        maximum NM tag [255].\n");
  fprintf(stderr, "     -v        verbose (print additional info for diagnosis).\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {

  int c;
  char *reffn = 0;
  char *reg = 0;
  char *infn = 0;
  char *outfn = 0;
  conf_t conf = {
    .step = 100000,
    .n_threads = 3,
    .min_base_qual = 20,
    .min_mapq = 40,
    .max_retention = 999999,
    .min_read_len = 10,
    .filter_qcfail = 1,
    .filter_secondary = 1,
    .filter_duplicate = 1,
    .filter_ppair = 1,
    .min_dist_end = 3,
    .max_nm = 255,
    .verbose = 0,
  };


  if (argc<2) return usage();
  while ((c=getopt(argc, argv, "i:o:r:g:q:e:b:t:n:m:l:cupvh"))>=0) {
    switch (c) {
    case 'i': infn = optarg; break;
    case 'r': reffn = optarg; break;
    case 'g': reg = optarg; break;
    case 'o': outfn = optarg; break;
    case 's': conf.step = atoi(optarg); break;
    case 'q': conf.n_threads = atoi(optarg); break;
    case 'b': conf.min_base_qual = atoi(optarg); break;
    case 't': conf.max_retention = atoi(optarg); break;
    case 'l': conf.min_read_len = atoi(optarg); break;
    case 'e': conf.min_dist_end = atoi(optarg); break;
    case 'c': conf.filter_secondary = 0; break;
    case 'u': conf.filter_duplicate = 0; break;
    case 'p': conf.filter_ppair = 0; break;
    case 'm': conf.min_mapq = atoi(optarg); break;
    case 'n': conf.max_nm = atoi(optarg); break;
    case 'v': conf.verbose = 1; break;
    case 'h': return usage();
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
  }

  if (!infn || !reffn) {
    usage();
    exit(1);
  }

  pileup_cytosine_main(reffn, infn, outfn, reg, &conf);
}
