/* convert bam to epiread format with supplied snp bed file
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
#include <zlib.h>
#include "bisc_utils.h"
#include "pileup.h"
#include "wzmisc.h"

typedef struct {
    bisc_common_t comm; /* common parameters across subcommands */
    bisc_threads_t bt;  /* multithreading parameters */
    meth_filter_t filt; /* methylation extraction filters */

    uint32_t epiread_reg_start;      /* first location of region provided to epiread */
    uint32_t epiread_reg_end;        /* final location of region provided to epiread */
    uint8_t  filter_empty_epiread:1; /* remove epireads that only have F/x/P */
    uint8_t  is_long_read:1;         /* data is from long read sequencing */
    uint8_t  epiread_old:1;          /* print old BISCUIT epiread format */
    uint8_t  epiread_pair:1;         /* pairwise output mode in epireads, doesn't mean "paired-end" */
    uint8_t  print_all_locations:1;  /* print all CpG and SNP locations in location column of epiread format */
} epiread_conf_t;

void epiread_conf_init(epiread_conf_t *conf);
