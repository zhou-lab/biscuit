/* extract barcodes from FASTQ file
 * 
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
 *
 */

#ifndef _BC_H_
#define _BC_H_

#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <zlib.h>

#include "kstring.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// Number of extra bytes to allocate
#define N_EXTRA 32

// uint8_t should be large enough for now, but if barcodes start to occur
// in locations further in from the start of the read or barcode lengths
// increase, then this will need to be bumped up in size
typedef struct {
    uint8_t  mate;      /* which read the barcode is on (1 or 2) */
    uint8_t  bc_start;  /* start position of barcode (1-based) */
    uint8_t  bc_length; /* length of barcode */
    char    *ofile;     /* prefix for output files, NULL writes to stdout */
} bc_conf_t;

static inline void bc_conf_init(bc_conf_t *conf) {
    conf->mate      = 1;
    conf->bc_start  = 1;
    conf->bc_length = 0;
    conf->ofile     = NULL;
}

gzFile setup_output(const char *ofile, uint8_t read_num);

void extract_barcodes(bc_conf_t *conf, kseq_t *ks1, kseq_t *ks2, gzFile oh1, gzFile oh2);

#endif /* _BC_H_ */
