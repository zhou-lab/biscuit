#include "indel_realign.h"


typedef struct insb_t insb_t;     /* insertion branch */
typedef struct delb_t delb_t;     /* deletion branch */

DEFINE_VECTOR(insb_v, insb_t)
DEFINE_VECTOR(delb_v, delb_t)

/* node in De Bruijn Graph */
typedef struct dbj_t {
  char c;
  insb_v *insb;
  delb_v *delb;
} dbj_t;

struct insb_t {
  char *s;
  int l;
  unsigned right;
};

struct delb_t {
  int l;                        /* deletion length */
  unsigned right;               /* right aligned pos */
}

#define bam_free1(b) free((b)->data)
KMEMPOOL_INIT(read, bam1_t, bam_free1)

typedef struct {
  dbj_graph_t *graph;
  samfile_t *sam;
  kmempool_t(read) *rmp;
  kmempool_t(dbj) *dbjmp;
  refseq_t *rs;
  bam1_t *b;
} detector_t;

void dbj_graph_put_insertion(uint32_t pos, char *seq, int l) {
  
}

void dbj_graph_put_deletion(uint32_t pos, unsigned l) {

  left = del_roll_left;
  right = del_roll_right
  
}

void log_indels(detector_t *dtt) {

  bam1_t *b = dtt->next;
  bam1_core_t *c = &b->core;
  uint32_t *cigar = bam1_cigar(b);

  uint32_t q=0, r=c->pos, op=0;
  uint32_t i, l;
  for (i=0; i<c->n_cigar; ++i) {
    op = bam_cigar_op(cigar[i]);
    l = bam_cigar_oplen(cigar[i]);

    switch (op) {
    case BAM_CMATCH:
      r += l;
      q += l;
      break;
    case BAM_CINS: {
      dbj_graph_put_insertion(r,l);
      q += l;
      break;
    }
    case BAM_CSOFT_CLIP:
      q += l;
      break;
    case BAM_CDEL: {
      dbj_graph_put_deletion(r,l);
      r += l;
      break;
    }
    case BAM_CREF_SKIP:
      r += l;
      break;
    default:
      fprintf(stderr, "Unknown cigar, %u\n", op);
      abort();
    }
  }
}

/* go through all reads mapped to first 1kb */
detector_t *init_detector(int tid, size_t pos, int pend, refseq_t *rs) {
  detector_t *dtt = (detector_t*) calloc(1,sizeof(detector_t));
  uint8_t regional;
  if (pend >= 0) {
    bam_iter_t iter = bam_iter_query(idx, tid, beg, pend);
    regional = 1;
  } else {
    bam_iter_t iter = bam_iter_query(idx, tid, beg, 1<<29);
    bam_iter_seek(sam->x.bam, iter);
    regional = 0;
  }

  dtt->rmp = kmp_init(read);
  dtt->rs = init_refseq(ref_fn, 10000, 10000);
  fetch_refseq(dtt->rs, target_name[tid], pos-1000, pos+2000);

  /* initialize local de-bruijn graph */
  dtt->graph = init_dbj_graph(3000);
  for (i=pos; i<pend; ++i) {
    dbj_t *n = ref_new_dbj_graph(dtt->graph, i);
    n->c = getbase_refseq(dtt->rs, i);
  }

  bam1_t *b = kmp_alloc(read, rmp);
  int ret;
  dtt->pos = pos;
  while (1) {
    if (regional) {
      if ((ret = bam_iter_read(sam->x.bam, iter, b))<0) {
        dtt->next = 0;
        break;
      }
    } else {
      if ((ret = bam_read1(sam->x.bam, b))<0) {
        dtt->next = 0;
        break;
      }
    }
    dtt->next = b;
    if (dtt->next->core.pos > dtt->pos || !dtt->next) break;
    log_indels(dtt);
  }
}
