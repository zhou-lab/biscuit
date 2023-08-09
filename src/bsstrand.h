/* Correct bisulfite strand information if it is very inconsistent with C2T/G2A count
 *
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
 *
 */

#ifndef _BSSTRAND_H_
#define _BSSTRAND_H_

#include <unistd.h>
#include "wzmisc.h"
#include "refcache.h"
#include "sam.h"
#include "bamfilter.h"

/* TAG_BSW      - bisulfite converted Watson
 * TAG_BSC      - bisulfite converted Crick
 * TAG_CONFLICT - conflicting evidences, coexistence of both C>T and G>A
 * TAG_UNKNOWN  - neither C>T or G>A exists */
typedef enum {TAG_BSW, TAG_BSC, TAG_CONFLICT, TAG_UNKNOWN} conversion_tag_t;

/* f - bisulfite converted Watson
 * r - bisulfite converted Crick
 * c - conflicting evidences, coexistence of both C>T and G>A
 * u - neither C>T or G>A exists */
static const char conversion_tags[4] = "frcu";
//const char conversion_tags[4] = "frcu";

typedef struct {
    uint8_t output_count;
    uint8_t correct_bsstrand;
} bsstrand_conf_t;

typedef struct {
    refcache_t *rs;
    int n_corr, n_mapped, n_unmapped;
    int confusion[16];
    int strandcnt[16];
    bsstrand_conf_t *conf;
} bsstrand_data_t;

/* f - bisulfite converted Watson
 * r - bisulfite converted Crick
 * c - conflicting evidences, coexistence of both C>T and G>A
 * u - neither C>T or G>A exists */
conversion_tag_t bam_tag_get_bsstrand(bam1_t *b);

int bsstrand_func(bam1_t *b, samFile *out, bam_hdr_t *header, void *data);

#endif /* _BSSTRAND_H_ */
