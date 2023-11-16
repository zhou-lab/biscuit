/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2023 Jacob.Morrison@vai.org
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

#include "wqueue.h"
#include "wvec.h"

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

#endif /* _BISC_UTILS_H_ */
