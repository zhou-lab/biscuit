/* utility functions for biscuit
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2023-2024 Jacob.Morrison@vai.org
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
**/

#ifndef _BISC_UTILS_H_
#define _BISC_UTILS_H_

#include <inttypes.h>
#include <ctype.h>

#include "wqueue.h"
#include "wvec.h"
#include "wzmisc.h"
#include "encode.h"

#include "hts.h"
#include "sam.h"

#include "refcache.h"

// Common parameters used in biscuit subcommands
typedef struct {
    uint8_t is_nome:1; /* input data is NOMe-seq */
    uint8_t verbose:1; /* print out extra information during processing */
} bisc_common_t;

static inline bisc_common_t bisc_common_init() {
    bisc_common_t out;

    out.is_nome = 0;
    out.verbose = 0;

    return out;
}

// Multithreading specific parameters
// TODO: can probably merge this into bisc_common_t
typedef struct {
    int step;      /* step size of window dispatching */
    int n_threads; /* number of processing threads */
} bisc_threads_t;

static inline bisc_threads_t bisc_threads_init() {
    bisc_threads_t out;

    out.step = 100000;
    out.n_threads = 3;

    return out;
}

// Parameters affecting methylation extraction
typedef struct {
    // Minimum values
    uint32_t min_base_qual;      /* minimum base quality */
    uint32_t min_read_len;       /* minimum read length */
    uint8_t  min_dist_end_5p;    /* minimum distance to 5' end of read to be included in methylation count */
    uint8_t  min_dist_end_3p;    /* minimum distance to 3' end of read to be included in methylation count */
    uint8_t  min_mapq;           /* minimum MAPQ score to include read */
    int      min_score;          /* minimum alignment score (from AS tag) */

    // Maximum values
    int      max_nm;             /* maximum number of non-cytosine-conversion mismatches (from NM tag) */
    uint32_t max_retention;      /* maximum cytosine retention in read */

    // Filter flags
    uint8_t  filter_ppair:1;     /* only use BAM_FPROPER_PAIR reads (hidden CLI option) */
    uint8_t  filter_secondary:1; /* don't include BAM_FSECONDARY reads (hidden CLI option) */
    uint8_t  filter_duplicate:1; /* don't include BAM_FDUP reads */
    uint8_t  filter_qcfail:1;    /* don't include BAM_FQCFAIL reads */
    uint8_t  filter_doublecnt:1; /* stop double-counting cytosine in overlapping mate reads */
} meth_filter_t;

static inline meth_filter_t meth_filter_init() {
    meth_filter_t out;

    out.min_base_qual = 20;
    out.min_read_len = 10;
    out.min_dist_end_5p = 3;
    out.min_dist_end_3p = 3;
    out.min_mapq = 40;
    out.min_score = 40;
    out.max_nm = 999999;
    out.max_retention = 999999;
    out.filter_ppair = 1;
    out.filter_secondary = 1;
    out.filter_duplicate = 1;
    out.filter_qcfail = 1;
    out.filter_doublecnt = 1;

    return out;
}

// Window blocks for processing regions
typedef struct {
    int64_t  block_id; /* ID number for window */
    int32_t  tid;      /* contig ID number of region */
    uint32_t beg;      /* beginning of region for window */
    uint32_t end;      /* end of region for window */
} window_t;

DEFINE_WQUEUE(window, window_t)

// Contig info
typedef struct {
    int32_t   tid;  /* contig ID number */
    char     *name; /* contig name */
    uint32_t  len;  /* contig length */
} target_t;

DEFINE_VECTOR(target_v, target_t);

static inline int compare_targets(const void *a, const void *b) {
    return strcmp(((target_t*)a)->name, ((target_t*)b)->name);
}

// Macro function definitions
#define max(a,b)                \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b);   \
     _a > _b ? _a : _b; })

#define min(a,b)                \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b);   \
     _a > _b ? _b : _a; })

#ifndef kroundup32
/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

// Functions related to bisulfite reads
uint32_t cnt_retention(refcache_t *rs, bam1_t *b, uint8_t bsstrand);
uint32_t get_mate_length(char *m_cigar);
uint8_t infer_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual);
uint8_t get_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual, int allow_u);

// Parse region string
static inline void biscuit_parse_region(const char *reg, void *hdr, int *tid, int *beg, int *end) {
    const char *q = hts_parse_reg(reg, beg, end);
    if (q) {
        char *tmp = (char*)malloc(q - reg + 1);
        strncpy(tmp, reg, q - reg);
        tmp[q - reg] = 0;
        *tid = bam_name2id(hdr, tmp);
        free(tmp);
    }
    else {
        // not parsable as a region, but possibly a sequence named "foo:a"
        *tid = bam_name2id(hdr, reg);
        *beg = 0; *end = INT_MAX;
    }
}

// Mutation and methylation codes
#define NSTATUS_METH 3
#define NSTATUS_BASE 7
extern const char nt256int8_to_methcode[NSTATUS_METH];
extern const char nt256int8_to_basecode[NSTATUS_BASE];
typedef enum {METH_RETENTION, METH_CONVERSION, METH_NA} status_meth_t;
typedef enum {BASE_A, BASE_C, BASE_G, BASE_T, BASE_N, BASE_Y, BASE_R} status_base_t;

// Cytosine context code (different representation of context when nome-seq)
#define NCONTXTS 6 /* not including NA */
typedef enum {CTXT_HCG, CTXT_HCHG, CTXT_HCHH, CTXT_GCG, CTXT_GCHG, CTXT_GCHH, CTXT_NA} cytosine_context_t;
extern const char *cytosine_context[NCONTXTS+1];
extern const char *cytosine_context_nome[NCONTXTS+1];

cytosine_context_t fivenuc_context(refcache_t *rs, uint32_t rpos, char rb, char *fivenuc);

// Record related structs and functions
#define RECORD_QUEUE_END -2
#define RECORD_SLOT_OBSOLETE -1

// TODO: figure out if there's a way to pull the pileup-only portions out
typedef struct {
    /* these parameters are used in pileup and epiread */
    int64_t   block_id; /* ID of block processed by thread */
    kstring_t s;        /* stores entries to print */
    int       tid;      /* contig ID number */

    /* these parameters are only used in pileup */
    double  *betasum_context; /* CG, CHG, CHH */
    int64_t *cnt_context;     /* number of Cs in each context */
} record_t;

DEFINE_VECTOR(record_v, record_t)
DEFINE_WQUEUE(record, record_t)

void pop_record_by_block_id(record_v *records, int64_t block_id, record_t *record);
void put_into_record_v(record_v *records, record_t rec);

// modBAM related functions
static inline float calculate_mod_probability(int qual) {
    if (qual < 0) {
        return -1.0;
    }
    return ((float)(qual) + 0.5) / 256.0;
}

static inline uint8_t is_modbam_cpg(uint16_t flag, int strand, int can_base, char qb, char rb, refcache_t *rs, uint32_t pos) {
    if (can_base == 'C' && strand == 0) {
        if (qb == 'G' && (flag & BAM_FREVERSE)) {
            if (rb == 'G' && pos-1 >= rs->beg && refcache_getbase_upcase(rs, pos-1) == 'C') {
                return 1;
            }
        } else if (qb == 'C' && !(flag & BAM_FREVERSE)) {
            if (rb == 'C' && pos+1 <= rs->end && refcache_getbase_upcase(rs, pos+1) == 'G') {
                return 1;
            }
        }
    } else if (can_base == 'G' && strand == 1) {
        if (qb == 'C' && (flag & BAM_FREVERSE)) {
            if (rb == 'C' && pos+1 <= rs->end && refcache_getbase_upcase(rs, pos+1) == 'G') {
                return 1;
            }
        } else if (qb == 'G' && !(flag & BAM_FREVERSE)) {
            if (rb == 'G' && pos-1 >= rs->beg && refcache_getbase_upcase(rs, pos-1) == 'C') {
                return 1;
            }
        }
    }

    return 0;
}

#endif /* _BISC_UTILS_H_ */
