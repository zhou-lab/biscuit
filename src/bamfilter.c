/** A simple filter of a bam file.
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021-2022 Jacob.Morrison@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rigbam
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

#include <stdlib.h>
#include "bamfilter.h"
#include "sam.h"
#include "wzmisc.h"

int bam_filter(char *ifn, char *ofn, char *reg, void *data, bam_filter_f func) {
  int ret = 0;
  samFile *in = sam_open(ifn, "rb");
  bam_hdr_t *header = sam_hdr_read(in);
  samFile *out = 0;
  if (ofn) {
    out = sam_open(ofn, "wb");
    if (!out) wzfatal("Cannot write bam %s.\n", ofn);
    if (sam_hdr_write(out, header) < 0) wzfatal("Cannot write bam header.\n");
  }

  if (reg) {

    hts_idx_t *idx = sam_index_load(in, ifn);
    if (idx == 0)
      wzfatal("[%s:%d] Random alignment retrieval only works for indexed BAM files.\n", __func__, __LINE__);

    int tid=0, beg, end;
    char *_reg = strdup(reg); /* a temporary string to modify */
    char *name_lim = (char *) hts_parse_reg(_reg, &beg, &end);
    if (name_lim) {
      char name_terminator = *name_lim;
      *name_lim = '\0';
      tid = bam_name2id(header, _reg);
      *name_lim = name_terminator;
    }
    free(_reg);
    if (tid < 0)
      wzfatal("[%s:%d] Region \"%s\" specifies an unknown reference name. \n", __func__, __LINE__, reg);

    bam1_t *b = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(idx, tid, beg, end);
    while ((ret = sam_itr_next(in, iter, b)) >= 0)
      func(b, out, header, data);
    hts_itr_destroy(iter);
    bam_destroy1(b);
    hts_idx_destroy(idx);

  } else {

    bam1_t *b = bam_init1();
    while ((ret = sam_read1(in, header, b)) >= 0)
      func(b, out, header, data);
    bam_destroy1(b);

  }

  if (out) sam_close(out);
  sam_close(in);
  bam_hdr_destroy(header);

  if (ret != -1) 				/* truncated is -2 */
    wzfatal("[%s:%d] Alignment retrieval failed due to truncated file\n", __func__, __LINE__);

  return ret;
}


