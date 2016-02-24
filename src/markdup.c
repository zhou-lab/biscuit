/********** mark duplicates ********
 * This utility requires bam file to be sorted before and after.
 *
 * reads are duplicate if they
 * 1. are mapped to the same coordinate
 * 2. have the same cigar string
 * 3. from the same bisulfite strand
 * 
 * From each insert group, one pick the read with the highest quality score.
 * 
 * run the script on a bam with duplicate marking will re-mark the duplicates
 * 
 * */

#include <stdlib.h>
#include "sam.h"
#include "khash.h"
#include "klist.h"
#include "wvec.h"
// #include "wstr.h"

void bam_sort_core_ext(int is_by_qname, const char *fn, const char *prefix, size_t _max_mem, int is_stdout, int n_threads, int level, int full_path);

typedef struct {
  uint32_t min_baseQ;		/* threshold for high quality base */
  uint8_t rmdup;
  uint8_t sort;
  uint32_t dup_cnt_pe;
  uint32_t cnt_pe;
  uint32_t dup_cnt_se;
  uint32_t cnt_se;
  uint32_t cnt_dangle;
  uint8_t verbose;
  uint8_t quiet;
  int max_isize;
  int mate_unmapped_as_se;
} mkconf_t;


typedef struct {
  bam1_t *b1;
  bam1_t *b2;
  uint32_t sum_qual;
  uint8_t is_dup:1;
  uint8_t is_PE:1;
  uint8_t bsstrand:1;
} __attribute__((__packed__)) insert_t;

DEFINE_VECTOR(insert_v, insert_t*)


int insert_hash_equal(insert_t *ins1, insert_t *ins2) {

  /* differentiate single-end from paired-end */
  if (ins1->is_PE!=ins2->is_PE) return 0;

  /* compare bsstrand */
  if (ins1->bsstrand != ins2->bsstrand) return 0;
  
  /* compare mate 1 */
  bam1_t *b11 = ins1->b1;
  bam1_t *b21 = ins2->b1;
  if (b11) {
    if (b21) {
      bam1_core_t *c11 = &b11->core;
      bam1_core_t *c21 = &b21->core;
      if (c11->tid != c21->tid) return 0;
      if (c11->pos != c21->pos) return 0;
    } else return 0;
  } else if (b21) return 0;

  /* compare mate 2 */
  bam1_t *b12 = ins1->b2;
  bam1_t *b22 = ins2->b2;
  if (b12) {
    if (b22) {
      bam1_core_t *c12 = &b12->core;
      bam1_core_t *c22 = &b22->core;
      if (c12->tid != c22->tid) return 0;
      if (c12->pos != c22->pos) return 0;
    } else return 0;
  } else if (b22) return 0;

  return 1;
}

#define insert_hash_func(insert) (khint32_t) ((insert)->b1 ? (insert)->b1->core.pos : (insert)->b2->core.pos)

KHASH_INIT(IGMap, insert_t*, insert_v*, 1, insert_hash_func, insert_hash_equal)

#define qname_hash_func(b) kh_str_hash_func(bam1_qname(b))

int read_match(bam1_t *b1, bam1_t *b2) {
  if ((strcmp(bam1_qname(b1), bam1_qname(b2)) == 0) &&
      (b1->core.mpos == b2->core.pos) &&
      (b1->core.pos == b2->core.mpos)) return 1;
  else return 0;
}

#define __mkdup_free_insert(key)	/* no-op */

KHASH_INIT(RIMap, bam1_t*, insert_t*, 1, qname_hash_func, read_match)
/* KHASH_MAP_INIT_STR(NameInsMap, insert_t*) */
#define bam_free1(b) free((b)->data)
KMEMPOOL_INIT(read, bam1_t, bam_free1)
KMEMPOOL_INIT(insert, insert_t, __mkdup_free_insert)

static inline int
sum_qual(const bam1_t *b)
{
  if (!b) return 0;
  int i, q;
  uint8_t *qual = bam1_qual(b);
  for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
  return q;
}

int cmp_insert_qual(const void *p1, const void *p2) {
  return (*((insert_t**)p2))->sum_qual - (*((insert_t**)p1))->sum_qual;
}

/* maximize the quality string of the
 * first insert by the quality string of the
 * second insert */
void maximize_qual1(insert_t *ins1, insert_t *ins2) {
  uint8_t *q1;
  uint8_t *q2;
  int x;
  if (ins1->b1) {
    q1 = bam1_qual(ins1->b1);
    q2 = bam1_qual(ins2->b1);
    for (x = 0; x<ins1->b1->core.l_qseq; ++x) if (q1[x] < q2[x]) q1[x] = q2[x];
  }
  if (ins1->b2) {
    q1 = bam1_qual(ins1->b2);
    q2 = bam1_qual(ins2->b2);
    for (x = 0; x<ins1->b2->core.l_qseq; ++x) if (q1[x] < q2[x]) q1[x] = q2[x];
  }
}

/* if all the high quality base of test
   agree with good, return 1
   else return 0 */
int
highqual_equal(insert_t *good, insert_t *test, uint32_t min_baseQ) {

  int32_t i;

  /* compare mate 1 */
  bam1_t *b11 = good->b1;
  bam1_t *b21 = test->b1;
  if (b11) {
    if (b21) {
      uint8_t *s11 = bam1_seq(b11);
      uint8_t *s21 = bam1_seq(b21);
      uint8_t *q21 = bam1_qual(b21);
      bam1_core_t *c11 = &b11->core;
      bam1_core_t *c21 = &b21->core;
      if (c11->l_qseq != c21->l_qseq) return 0;
      for (i=0; i<c21->l_qseq; ++i) {
        if (bam1_seqi(s11, i) != bam1_seqi(s21, i) &&
            q21[i] > min_baseQ)
          return 0;
      }
    } else return 0;
  } else if (b21) return 0;

  /* compare mate 2 */
  bam1_t *b12 = good->b2;
  bam1_t *b22 = test->b2;
  if (b12) {
    if (b22) {
      uint8_t *s12 = bam1_seq(b12);
      uint8_t *s22 = bam1_seq(b22);
      uint8_t *q22 = bam1_qual(b22);
      bam1_core_t *c12 = &b12->core;
      bam1_core_t *c22 = &b22->core;
      if (c12->l_qseq != c22->l_qseq) return 0;
      for (i=0; i<c22->l_qseq; ++i) {
        if (bam1_seqi(s12, i) != bam1_seqi(s22, i) &&
            q22[i] > min_baseQ)
          return 0;
      }
    } else return 0;
  } else if (b22) return 0;

  return 1;
  
}

void resolve_dup(khash_t(IGMap) *igm, samfile_t *out,
                 kmempool_t(read) *rmp,
                 kmempool_t(insert) *imp,
                 mkconf_t *conf) {

  uint8_t rmdup = conf->rmdup;
  khint_t k;
  for (k=kh_begin(igm); k<kh_end(igm); ++k) {
    if (kh_exist(igm, k)) {
      insert_v *ig = kh_val(igm, k);

      /* keep the read with the best sum of quality, 
       * and mark rest as duplicates */
      unsigned i;
      unsigned bestqual;
      insert_t *bestins=0;
      for (i=0; i<ig->size; ++i) {
        insert_t *ins = get_insert_v(ig, i);
        ins->sum_qual = sum_qual(ins->b1) + sum_qual(ins->b2);
        if (i==0 || ins->sum_qual>bestqual) {
          bestins = ins;
          bestqual=ins->sum_qual;
        }
      }

      bestins->is_dup = 0;
      for (i=0; i<ig->size; ++i) {
        insert_t *ins = get_insert_v(ig, i);
        if (ins != bestins) ins->is_dup = 1;
      }

      /* dump the insert group */
      if (conf->verbose > 5) {
        if (ig->size) printf("Insert group: %s\n", ig->buffer[0]->b1?bam1_qname(ig->buffer[0]->b1):bam1_qname(ig->buffer[0]->b2));
      }
      for (i=0; i<ig->size; ++i) {
        insert_t *ins = get_insert_v(ig, i);
        if (ins->is_dup) {
          if (ins->is_PE) conf->dup_cnt_pe += 2;
          else ++conf->dup_cnt_se;
          if (ins->b1) ins->b1->core.flag |= BAM_FDUP;
          if (ins->b2) ins->b2->core.flag |= BAM_FDUP;
        }
        if (!(ins->is_dup && rmdup)) {
          if (ins->b1) {
            samwrite(out, ins->b1);
            kmp_free(read, rmp, ins->b1);
          }
          if (ins->b2) {
            samwrite(out, ins->b2);
            kmp_free(read, rmp, ins->b2);
          }
        }
        kmp_free(insert, imp, ins);
      }
      free_insert_v(ig);
    }
  }
  kh_clear(IGMap, igm);
}

static int8_t bam_get_bsstrand(bam1_t *b) {
  uint8_t *s;

  s = bam_aux_get(b, "ZS");     /* bsmap flag */
  if (s) {
    s++;
    if (*s == '+') return 0;
    else if (*s == '-') return 1;
    else {
      fprintf(stderr, "[%s:%d] Unknown ZS strand tag: %c\n", __func__, __LINE__, *s);
      fflush(stderr);
      exit(1);
    }
  }

  s = bam_aux_get(b, "YD");     /* bwa-meth flag */
  if (s) {
    s++;
    if (*s == 'f') return 0;
    else if (*s == 'r') return 1;
    else if (*s == 'u') return -1;
  }

  s = bam_aux_get(b, "XG");     /* bismark flag */
  if (s) {
    s++;
    if (strcmp((char*)s, "CT")==0) return 0;
    else if (strcmp((char*)s, "GA")) return 1;
  }

  /* if no available flag information */
  fprintf(stderr, "[%s:%d] Warning: bisulfite strand information missing, treat as BSW: %s\n", __func__, __LINE__, bam1_qname(b));
  fflush(stderr);
  return 0;
}

static void flush_dangling_reads(khash_t(RIMap) *rim, kmempool_t(read) *rmp, kmempool_t(insert) *imp, samfile_t *out, mkconf_t *conf) {
  khint_t k;
  for (k=kh_begin(rim); k<kh_end(rim); ++k) {
    if (kh_exist(rim, k)) {
      insert_t *ins = kh_val(rim, k);
      if (ins->b1) {
        if (conf->verbose > 5) {
          fprintf(stderr, "[%s:%d] Warning: found dangling reads: %s\n", __func__, __LINE__, bam1_qname(ins->b1));
          fflush(stderr);
        }
        samwrite(out, ins->b1);
        kmp_free(read, rmp, ins->b1);
        conf->cnt_dangle++;
      }
      if (ins->b2) {
        if (conf->verbose > 5) {
          fprintf(stderr, "[%s:%d] Warning: found dangling reads: %s\n", __func__, __LINE__, bam1_qname(ins->b2));
          fflush(stderr);
        }
        samwrite(out, ins->b2);
        kmp_free(read, rmp, ins->b2);
        conf->cnt_dangle++;
      }
      kmp_free(insert, imp, ins);
      kh_del(RIMap, rim, k);
    }
  }
}

int mark_dup(char *bam_in_fn, char *bam_out_fn, mkconf_t *conf) {

  char *bam_out0_fn;
  if (conf->sort) {
    bam_out0_fn = malloc(strlen(bam_out_fn)+5);
    sprintf(bam_out0_fn, "%s.tmp", bam_out_fn);
  } else {
    bam_out0_fn = bam_out_fn;
  }

  samfile_t *in = samopen(bam_in_fn, "rb", 0);
  samfile_t *out = samopen(bam_out0_fn, "wb", in->header);

  int last_tid = -1, last_pos = -1;

  khash_t(IGMap) *igm = kh_init(IGMap);
  khash_t(RIMap) *rim = kh_init(RIMap);

  kmempool_t(read) *rmp = kmp_init(read);
  kmempool_t(insert) *imp = kmp_init(insert);

  unsigned cnt = 0;
  bam1_t *b = kmp_alloc(read, rmp);
  while (samread(in, b)>=0) {
    cnt++;
    /* fprintf(stderr, "mtid: %d isize: %d\n", b->core.mtid, b->core.isize); */
    if (!conf->quiet && (cnt & 0xFFF)==0) {
      fprintf(stderr, "\r[%s] parsed %u reads.", __func__, cnt);
      fflush(stderr);
    }

    bam1_core_t *c = &b->core;
    c->flag &= ~BAM_FDUP;       /* remove existing duplication info */

    /* process insert group */
    if (c->tid != last_tid || c->pos != last_pos) {
      resolve_dup(igm, out, rmp, imp, conf);

      if (c->tid != last_tid) {
        flush_dangling_reads(rim, rmp, imp, out, conf);
      }

      last_tid = c->tid; last_pos = c->pos;
    }

    /* skip unmapped and secondary */
    if (c->flag & BAM_FUNMAP || c->flag & BAM_FSECONDARY) {
      samwrite(out, b);
      continue;
    }

    int ret;
    /* if mate is unmapped in PE, treat as SE */
    if ((c->flag & BAM_FPAIRED) && !((c->flag&BAM_FMUNMAP) && conf->mate_unmapped_as_se)) { /* PE */

      ++conf->cnt_pe;
      khint_t e = kh_put(RIMap, rim, b, &ret);
      insert_t *ins;
      if (ret) {			/* empty */
        ins = kmp_alloc(insert, imp);
        memset(ins, 0, sizeof(insert_t));
        kh_val(rim, e) = ins;
      } else {
        ins = kh_val(rim, e);
      }

      ins->is_PE = 1;
      if (c->flag & BAM_FREAD1) {
        ins->b1 = b;
        b = kmp_alloc(read, rmp);
      } else if (c->flag & BAM_FREAD2) {
        ins->b2 = b;
        b = kmp_alloc(read, rmp);
      } else {
        fprintf(stderr, "[%s:%d] read position undecided, skip.\n", __func__, __LINE__);
      }

      /* Once the insert is matched or mate-unmapped or mate problematic,
       * push into insert group.
       * Dangling reads are handled at the end. */
      if ((ins->b1 && ins->b2)
          || (ins->b1 && ins->b1->core.tid != ins->b1->core.mtid) || (ins->b2 && ins->b2->core.tid != ins->b2->core.mtid)
          || (ins->b1 && abs(ins->b1->core.isize) >conf->max_isize) || (ins->b2 && abs(ins->b2->core.isize) >conf->max_isize)
          || (ins->b1 && (ins->b1->core.flag&BAM_FMUNMAP)) || (ins->b2 && (ins->b2->core.flag&BAM_FMUNMAP))) {

        /* find bsstrand of the whole insert */
        int bs1=ins->b1 ? bam_get_bsstrand(ins->b1) : -1;
        int bs2=ins->b2 ? bam_get_bsstrand(ins->b2) : -1;
        if (bs1 == bs2 && bs1 >=0) {
          ins->bsstrand = (unsigned) bs1;
        } else {
          /* fix mate's bsstrand */
          if (bs1 < 0 && bs2 >= 0) ins->bsstrand = (unsigned) bs2;
          else if (bs2 < 0 && bs1 >= 0) ins->bsstrand = (unsigned) bs1;
          else if (bs1 < 0 && bs2 < 0) {
            /* in rare circumstances, one is unmapped, the other has
             * undetermined bsstrand. then assume T-rich conversion strand */
            if (conf->verbose>0) {
              fprintf(stderr, "[%s:%d] No valid BS strand info: %s\n", __func__, __LINE__, bam1_qname(ins->b1));
              fflush(stderr);
            }
            ins->bsstrand = 0;
          } else {
            if (conf->verbose>0) {
              fprintf(stderr, "[%s:%d] Warning: inconsistent bisulfite strand between mate reads of %s\n", __func__, __LINE__, bam1_qname(ins->b1));
              fflush(stderr);
            }
            ins->bsstrand = (unsigned) (sum_qual(ins->b1)>sum_qual(ins->b2) ? bs1 : bs2);
          }
        }

        insert_v *ig;
        khint_t k = kh_put(IGMap, igm, ins, &ret);

        if (ret) {		/* empty */
          ig = init_insert_v(2);
          kh_val(igm, k) = ig;
        } else {
          ig = kh_val(igm, k);
        }
        push_insert_v(ig, ins);
        /* remove matched reads from read-insert map */
        kh_del(RIMap, rim, e);
      }

    } else {			 /* SE */

      ++conf->cnt_se;
      insert_t *ins = kmp_alloc(insert, imp);
      memset(ins, 0, sizeof(insert_t));
      ins->is_PE = 0;
      ins->b1 = b;
      int bs = bam_get_bsstrand(ins->b1);
      if (bs <= 0) ins->bsstrand = 0;
      else ins->bsstrand = 1;
      b = kmp_alloc(read, rmp);

      insert_v *ig;
      khint_t k = kh_put(IGMap, igm, ins, &ret);
      if (ret) {		/* empty */
        ig = init_insert_v(2);
        kh_val(igm, k) = ig;
      } else {
        ig = kh_val(igm, k);
      }
      push_insert_v(ig, ins);
    }

  }
  resolve_dup(igm, out, rmp, imp, conf);
  flush_dangling_reads(rim, rmp, imp, out, conf);

  fprintf(stderr, "\r[%s] parsed %u reads\n", __func__, cnt);
  fprintf(stderr, "[%s] marked %d duplicates from %d paired-end reads (%.3g%%)\n", __func__, conf->dup_cnt_pe, conf->cnt_pe, (double) (conf->dup_cnt_pe) / conf->cnt_pe * 100);
  fprintf(stderr, "[%s] marked %d duplicates from %d single-end reads (%.3g%%)\n", __func__, conf->dup_cnt_se, conf->cnt_se, (double) (conf->dup_cnt_se) / conf->cnt_se * 100);
  fprintf(stderr, "[%s] identified %d dangling paired-end reads (%.3g%%)\n", __func__, conf->cnt_dangle, (double) conf->cnt_dangle / cnt * 100);
  fflush(stderr);
  
  samclose(in);
  samclose(out);

  kmp_free(read, rmp, b);
  kmp_destroy(read, rmp);
  kmp_destroy(insert, imp);
  kh_destroy(IGMap, igm);
  kh_destroy(RIMap, rim);

  if (conf->sort) {
    fprintf(stderr, "[%s] sorting after mkdup\n", __func__);
    bam_sort_core_ext(0, bam_out0_fn, bam_out_fn, 768<<20, 0, 0, -1, 1);
    bam_index_build(bam_out_fn);
    remove(bam_out0_fn);
    free(bam_out0_fn);
  }

  return 0;
}

/* without sorting */
int mark_dup_nosort(char *bam_in_fn, char *bam_out_fn) {
  mkconf_t conf = {
    .min_baseQ = 30,
    .rmdup = 0,
    .sort = 0,
    .dup_cnt_se = 0,
    .dup_cnt_pe = 0,
    .cnt_se = 0,
    .cnt_pe = 0,
    .cnt_dangle = 0,
    .verbose = 0,
    .quiet = 0,
    .max_isize = 10000,
    .mate_unmapped_as_se = 0,
  };

  return mark_dup(bam_in_fn, bam_out_fn, &conf);
}

static int usage(mkconf_t *conf) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: biscuit markdup [options] <in.bam> <out.bam>\n");
  fprintf(stderr, "       <in.bam> must be sorted and indexed.\n");
  fprintf(stderr, "       <out.bam> is default to <in.bam.mkdup> if not provided.\n");
  fprintf(stderr, "Input options:\n");
  fprintf(stderr, "     -l INT    maximum insert size for paired end duplicate-marking [%d]\n", conf->max_isize);
  fprintf(stderr, "     -b INT    minimum base quality [30].\n");
  fprintf(stderr, "     -r        toggle to remove marked duplicate.\n");
  fprintf(stderr, "     -s        toggle to turn OFF sorting and indexing after marking duplicates.\n");
  fprintf(stderr, "     -u        treat mate-unmapped paired-end reads as single-end reads\n");
  fprintf(stderr, "     -q        quiet (log friendly)\n");
  fprintf(stderr, "     -v        verbose level [%d].\n", conf->verbose);
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main_markdup(int argc, char *argv[]) {

  mkconf_t conf = {
    .min_baseQ = 30,
    .rmdup = 0,
    .sort = 1,
    .dup_cnt_se = 0,
    .dup_cnt_pe = 0,
    .cnt_se = 0,
    .cnt_pe = 0,
    .verbose = 0,
    .quiet = 0,
    .max_isize = 10000,
    .mate_unmapped_as_se = 0,
  };

  int c;
  if (argc < 2) return usage(&conf);
  while ((c = getopt(argc, argv, "b:l:v:rsuqh")) >= 0) {
    switch (c) {
    case 'b': conf.min_baseQ = atoi(optarg); break;
    case 'l': conf.max_isize = atoi(optarg); break;
    case 'v': conf.verbose = atoi(optarg); break;
    case 'r': conf.rmdup = 1; break;
    case 's': conf.sort = 0; break;
    case 'u': conf.mate_unmapped_as_se = 0; break;
    case 'q': conf.quiet = 1; break;
    case 'h': return usage(&conf);
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "[%s:%d] missing input.\n", __func__, __LINE__);
    exit(1);
  }

  char *in_fn=0, *out_fn=0;
  in_fn = argv[optind];
  if (optind+1 < argc) out_fn = strdup(argv[optind+1]);
  else {
    out_fn = calloc(strlen(in_fn)+7, 1);
    strcpy(out_fn, in_fn);
    strcat(out_fn, ".mkdup");
  }

  mark_dup(in_fn, out_fn, &conf);

  free(out_fn);

  return 0;
}
