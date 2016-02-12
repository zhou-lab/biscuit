
#include <ctype.h>
#include <stdlib.h>
#include <inttypes.h>
#include <libgen.h>
#include "wqueue.h"
#include "encode.h"
#include "sam.h"
#include "refseq.h"
#include "kstring.h"
#include "wvec.h"
#include "stats.h"
#include "biscuit.h"

#define bscall(b, pos) bam_nt16_rev_table[bam1_seqi(bam1_seq(b), pos)]

typedef struct {
  int step;
  int n_threads;
  int min_cov;
  int bsrate_max_pos;
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
  uint8_t noheader:1;
  double error;
  double mu;
  double contam;
  double prior0;
  double prior1;
  double prior2;
  uint8_t verbose;
  int is_nome;
} conf_t;

void conf_init(conf_t *conf);

typedef struct {
  int64_t block_id;
  int32_t tid;
  uint32_t beg, end;
} window_t;

DEFINE_WQUEUE(window, window_t)

/* mutation-methylation code */
extern const char nt256int8_to_mutcode[6];
typedef enum {BSS_MA, BSS_MC, BSS_MG, BSS_MT,
              BSS_MY, BSS_MR, BSS_RETENTION, BSS_CONVERSION, BSS_N} status_t;

/* cytosine context */
#define NCONTXTS 6          /* not including NA */
typedef enum {CTXT_HCG, CTXT_HCHG, CTXT_HCHH,
              CTXT_GCG, CTXT_GCHG, CTXT_GCHH,
              CTXT_NA} cytosine_context_t;
extern const char *cytosine_context[];

typedef struct {
  uint8_t sid;			/* which sample */
  uint8_t bsstrand:1;
  uint8_t qual:7;
  uint8_t strand:1;
  uint16_t qpos;
  uint8_t cnt_ret;
  uint16_t rlen;                /* read length */
  char qb;
  status_t stat;                /* code from mut-met status table */
} __attribute__((__packed__)) pileup_data_t;

DEFINE_VECTOR(pileup_data_v, pileup_data_t)

#define RECORD_QUEUE_END -2
#define RECORD_SLOT_OBSOLETE -1

typedef struct bsrate_t {
  int m;
  int *ct_unconv, *ct_conv, *ga_unconv, *ga_conv;
  int *ct_unconv_m, *ct_conv_m, *ga_unconv_m, *ga_conv_m;
} bsrate_t;

static inline void bsrate_init(bsrate_t *b, int m) {
  b->m = m;
  b->ct_unconv = calloc(m, sizeof(int));
  b->ct_conv = calloc(m, sizeof(int));
  b->ga_unconv = calloc(m, sizeof(int));
  b->ga_conv = calloc(m, sizeof(int));

  /* conversion based on chrM */
  b->ct_unconv_m = calloc(m, sizeof(int));
  b->ct_conv_m = calloc(m, sizeof(int));
  b->ga_unconv_m = calloc(m, sizeof(int));
  b->ga_conv_m = calloc(m, sizeof(int));
}

static inline void bsrate_free(bsrate_t *b, int n_bams) {

  int sid;
  for (sid = 0; sid < n_bams; ++sid) {
    free(b[sid].ct_unconv);
    free(b[sid].ct_conv);
    free(b[sid].ga_unconv);
    free(b[sid].ga_conv);
    free(b[sid].ct_unconv_m);
    free(b[sid].ct_conv_m);
    free(b[sid].ga_unconv_m);
    free(b[sid].ga_conv_m);
  }
}

typedef struct {
  int64_t block_id;
  kstring_t s;                  /* vcf record */

  /* coverage */
  int tid;
  int64_t l;
  int64_t *n, *n_uniq;             /* length, base coverage, unique base coverage */

  /* methlevelaverages, [beta sum, cnt] */
  /* dim = NCONTXTS * n_bams */
  double *betasum_context;       /* CG, CHG, CHH */
  int64_t *cnt_context;
  
  /* bsrate */
  bsrate_t *b;
} record_t;

DEFINE_VECTOR(record_v, record_t)

DEFINE_WQUEUE(record, record_t)

void pop_record_by_block_id(record_v *records, int64_t block_id, record_t *record);
void put_into_record_v(record_v *records, record_t rec);

uint32_t cnt_retention(refseq_t *rs, bam1_t *b, uint8_t bsstrand);

uint8_t infer_bsstrand(refseq_t *rs, bam1_t *b, uint32_t min_base_qual);

uint8_t get_bsstrand(refseq_t *rs, bam1_t *b, uint32_t min_base_qual);

typedef struct {
  int32_t tid;
  char *name;
  uint32_t len;
} target_t;

DEFINE_VECTOR(target_v, target_t);

static inline int compare_targets(const void *a, const void *b) {
  return strcmp(((target_t*)a)->name, ((target_t*)b)->name);
}

#define mutcode(a) (nt256char_to_nt256int8_table[(uint8_t)a])

cytosine_context_t fivenuc_context(refseq_t *rs, uint32_t rpos, char rb, char *fivenuc);

#define min(a,b) ((a)<(b)?(a):(b))

void pileup_genotype(int cref, int altsupp, conf_t *conf, char gt[4], double *_gl0, double *_gl1, double *_gl2, double *_gq);
int reference_supp(int cnts[9]);
void allele_supp(char rb, int cref, int cm1, int cm2, int cnts[9], kstring_t *s);

static inline int compare_supp(const void *a, const void *b)
{
  return ((*(uint32_t*)b)>>4) - ((*(uint32_t*)a)>>4);
}

typedef struct {
  wqueue_t(record) *q;
  int n_bams;
  char **bam_fns;
  char *outfn;
  char *statsfn;
  char *header;
  target_v *targets;
  conf_t *conf;
} writer_conf_t;

void *write_func(void *data);


#ifndef kroundup32
/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif
