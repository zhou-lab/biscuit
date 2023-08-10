/* generate QC files
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2021-2023 Jacob.Morrison@vai.org
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

#ifndef _QC_H_
#define _QC_H_

#include <errno.h>
#include "wstr.h"
#include "bamfilter.h"
#include "bsstrand.h"
#include "bsconv.h"
#include "cinread.h"

// Variable Definitions
// -----------------------------------------------------------------------------
// MAPQ values run from 0 to 60, inclusive
// Set array limit at 62 to include unmapped (index 61) and MAPQs from 0 to 60
#define N_MAPQ 61

// There can be a large range of insert sizes
// Limit to insert sizes <= ISIZE
#define ISIZE 1000

// Define optimal MAPQ minimum
#define OPTIMAL_MAPQ 40

// Define output file tag lengths
#define LEN_SAMP 256
#define LEN_MAPQ 15
#define LEN_ISIZE 16
#define LEN_DUP 15
#define LEN_STRAND 17
#define LEN_BSCONV 28
#define LEN_CGREADPOS 26
#define LEN_CHREADPOS 26

// Type Definitions
// -----------------------------------------------------------------------------
typedef struct {
    FILE *f_out_mapq;
    FILE *f_out_isize;
    FILE *f_out_dup;
    FILE *f_out_strand;
    FILE *f_out_bsconv;
    FILE *f_out_cgreadpos;
    FILE *f_out_chreadpos;
} qc_outfiles_t;

typedef struct {
    int single_end;
    qc_outfiles_t out;
} qc_conf_t;


// Function Definitions
// -----------------------------------------------------------------------------
int  process_qc(char *input_bam, qc_conf_t conf, bsstrand_data_t *data_bsstrand, bsconv_data_t *data_bsconv, cinread_data_t *data_cinread_cg, cinread_data_t *data_cinread_ch);
void format_mapq_table(FILE *fname, unsigned int *values);
void format_isize_table(FILE *fname, unsigned int *values, unsigned int count);
void format_dup_report(FILE *fname, unsigned int all_tot, unsigned int all_dup, unsigned int q40_tot, unsigned int q40_dup);
void format_strand_report(FILE *fname, bsstrand_data_t *data);
void format_bsconv_report(FILE *fname, bsconv_data_t *data);
void format_readpos_report(FILE *fname, cinread_data_t *data, char *type);

static void close_outfiles(qc_outfiles_t out, qc_conf_t conf) {
    fclose(out.f_out_mapq);
    fclose(out.f_out_dup);
    fclose(out.f_out_strand);
    fclose(out.f_out_bsconv);
    fclose(out.f_out_cgreadpos);
    fclose(out.f_out_chreadpos);
    if (conf.single_end == 0)
        fclose(out.f_out_isize);
}

static void create_filename(char *out, char *samp, char *tag) {
    strcpy(out, samp);
    strcat(out, tag);
}

#endif /* _QC_H_ */
