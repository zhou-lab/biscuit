/* Print C in reads as long form
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021      Jacob.Morrison@vai.org
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
 */ 

#include <unistd.h>
#include <errno.h>
#include "wstr.h"
#include "wzmisc.h"
#include "refcache.h"
#include "sam.h"
#include "bamfilter.h"
#include "pileup.h"

static const char *tp_names[] = {
    "QNAME",     // read name
    "QPAIR",     // which read in pair
    "STRAND",    // forward or reverse strand
    "BSSTRAND",  // which original strand the read derives from
    "MAPQ",      // MAPQ score
    "QBEG",      // read start position
    "QEND",      // read end position
    "CHRM",      // chromosome
    "CRPOS",     // cytosine position on reference
    "CGRPOS",    // CpG position on reference (-1 if not applicable)
    "CQPOS",     // cytosine position on read
    "CRBASE",    // cytosine reference base
    "CCTXT",     // cytosine context, strand flipped
    "CQBASE",    // base called on read
    "CRETENTION" // retention (R) or conversion (C)
};

typedef enum {
    TP_QNAME,     // read name
    TP_QPAIR,     // which read in pair
    TP_STRAND,    // forward (+) or reverse (-) strand
    TP_BSSTRAND,  // which original strand the read derives from (+, BSW) or (-, BSC)
    TP_MAPQ,      // MAPQ score
    TP_QBEG,      // read start position
    TP_QEND,      // read end position
    TP_CHRM,      // chromosome
    TP_CRPOS,     // cytosine position on reference
    TP_CGRPOS,    // CpG position on reference (-1 if not applicable)
    TP_CQPOS,     // cytosine position on query
    TP_CRBASE,    // cytosine reference
    TP_CCTXT,     // cytosine context, strand flipped
    TP_CQBASE,    // base called on read
    TP_CRETENTION // retention (R) or conversion (C)
} __tp_name_t;

static const char *tgt_names[] = {"c", "cg", "ch", "hcg", "gch", "hch"};

typedef enum {SL_C, SL_CG, SL_CH, SL_HCG, SL_GCH, SL_HCH} __tgt_name_t; // SL_ select target

typedef struct {
    int n_tp_names;
    __tp_name_t *tp_names;
    __tgt_name_t tgt;
    FILE *out;
    int skip_secondary;
    int skip_printing;
} cinread_conf_t;

#define CIN_N_READS 2
#define CIN_READ_LEN 301
#define CIN_N_RET_STATES 3

typedef struct cinread_data_t {
    refcache_t *rs;
    FILE *output;
    cinread_conf_t *conf;
    int counts[CIN_N_READS][CIN_READ_LEN][CIN_N_RET_STATES];
} cinread_data_t;

int cinread_func(bam1_t *b, samFile *out, bam_hdr_t *hdr, void *data);
