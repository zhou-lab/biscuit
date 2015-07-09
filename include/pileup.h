
#include <ctype.h>
#include <stdlib.h>
#include "wqueue.h"
#include "encode.h"
#include "sam.h"
#include "refseq.h"
#include "kstring.h"
#include "wvec.h"
#include "stats.h"

#define bscall(b, pos) bam_nt16_rev_table[bam1_seqi(bam1_seq(b), pos)]

typedef struct {
  int64_t block_id;
  int32_t tid;
  uint32_t beg, end;
} window_t;

DEFINE_WQUEUE(window, window_t)

extern const char nt256int8_to_mutcode[6];
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

#define RECORD_QUEUE_END -2
#define RECORD_SLOT_OBSOLETE -1

typedef struct {
  int64_t block_id;
  kstring_t s;
} record_t;

DEFINE_VECTOR(record_v, record_t)

DEFINE_WQUEUE(record, record_t)


void remove_record_by_block_id(record_v *records, int64_t block_id, record_t *record);
void put_into_record_v(record_v *records, record_t rec);

uint32_t cnt_retention(refseq_t *rs, bam1_t *b, uint8_t bsstrand);

typedef struct {
  wqueue_t(record) *q;
  char *outfn;
} writer_conf_t;

void *write_func(void *data);

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

void fivenuc_context(refseq_t *rs, uint32_t rpos, kstring_t *s, char rb);

