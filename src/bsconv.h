/* annotate bisulfite conversion of reads
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

#ifndef _BSCONV_H_
#define _BSCONV_H_

#include <unistd.h>
#include <errno.h>
#include "wzmisc.h"
#include "refcache.h"
#include "sam.h"
#include "bamfilter.h"
#include "pileup.h"
#include "encode.h"

typedef struct {
    int max_cph;
    float max_cph_frac;
    float max_cpy_frac;
    int max_cpa;
    int max_cpc;
    int max_cpt;
    int filter_u;
    int show_filtered;
    int print_in_tab;
    int no_printing;
} bsconv_conf_t;

typedef struct bsconv_data_t {
    refcache_t *rs;
    bsconv_conf_t *conf;
    int n;
    int n_filtered;
    uint64_t retn_conv_counts[9];
} bsconv_data_t;

int bsconv_func(bam1_t *b, samFile *out, bam_hdr_t *hdr, void *data);

#endif /* _BSCONV_H_ */
