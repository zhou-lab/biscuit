/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021-2023 Jacob.Morrison@vai.org
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

#ifndef _PILEUP_H_
#define _PILEUP_H_

#include <unistd.h>
#include <ctype.h>
#include <stdlib.h>
#include <inttypes.h>
#include <libgen.h>

#include "wqueue.h"
#include "wzmisc.h"
#include "encode.h"
#include "wvec.h"
#include "stats.h"

#include "sam.h"
#include "hts.h"
#include "kstring.h"

#include "refcache.h"
#include "biscuit.h"
#include "bisc_utils.h"

typedef struct {
    bisc_common_t comm; /* common parameters across subcommands */
    bisc_threads_t bt;  /* multithreading parameters */
    meth_filter_t filt; /* methylation extraction filters */

    uint8_t ambi_redist:1; /* redistribute ambiguous (R/Y) calls in SNP genotyping */
    uint8_t somatic:1;     /* call somatic mutation */
    double error;          /* error rate */
    double mu;             /* mutation rate */
    double mu_somatic;     /* somatic mutation rate*/
    double contam;         /* contamination rate */
    double prior0;         /* prior probability for no variant (defined as 1 - prior1 - prior2) */
    double prior1;         /* prior probability for heterozygous variant */
    double prior2;         /* prior probability for homozygous variant */
} pileup_conf_t;

void pileup_conf_init(pileup_conf_t *conf);

typedef struct {
    uint8_t sid;        /* which sample */
    uint8_t bsstrand:1; /* bisulfite strand */
    uint8_t qual:7;     /* base quality */
    uint8_t strand:1;   /* read stand */
    uint16_t qpos;      /* position on read */
    uint8_t cnt_ret;    /* count of retention of entire read */
    uint16_t rlen;      /* read length */
    char qb;            /* query base */
    uint8_t stat;       /* (status_base_t << 4) | status_meth_t */
} __attribute__((__packed__)) pileup_data_t;

DEFINE_VECTOR(pileup_data_v, pileup_data_t)

void pileup_genotype(int cref, int altsupp, pileup_conf_t *conf, char gt[4], double *_gl0, double *_gl1, double *_gl2, double *_gq);

static inline int compare_supp(const void *a, const void *b) {
    return ((*(uint32_t*)b)>>4) - ((*(uint32_t*)a)>>4);
}

void *write_func(void *data);

#endif /* _PILEUP_H_ */
