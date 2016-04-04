/**
 * copy number segmentation
 * 
 * The MIT License (MIT)
 * Copyright (c) 2016 Wanding.Zhou@vai.org
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

#include "wzcbs.h"
#include "wzbed.h"



int main_copyseg(int argc, char *argv[])
{
  int c; int whole_ref = 0;
  char *bed_fn = NULL; char *region;
  int wsize = 500;
  while ((c = getopt(argc, argv, "g:r:wV:h")) >= 0) {
    switch (c) {
    case 'g': region = atoi(optarg); break;
    case 'r': bed_fn = atoi(optarg); break;
    case 'w': whole_ref = 1; break;
    case 'V': conf.verbose = atoi(optarg); break;
    case 'h': {
      fprintf(stderr, "\n");
      fprintf(stderr, "Usage: copyseg [options] input.bam \n");
      fprintf(stderr, "Input options:\n");
      fprintf(stderr, "     -r        bed file of targets [none]\n");
      fprintf(stderr, "     -g        target region [none]\n");
      fprintf(stderr, "     -s        step size [%d]\n", wsize);
      fprintf(stderr, "     -w        whole reference genome as given in bam header\n");
      fprintf(stderr, "     -V INT    verbose level [%d].\n", conf.verbose);
      fprintf(stderr, "     -h        this help.\n");
      fprintf(stderr, "\n");
      return 1;
    }
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      fflush(stderr);
      exit(1);
      break;
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "[%s:%d] Please provide input bam. Abort.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  samfile_t *in = samopen(argv[optind], "rb", 0);
  target_v *targets = bamheader2targets(in->header);

  if (bed_fn) {
    bed1_v *beds = bed_read(bed_fn, targets);
  } else if (region) {
    bed1_v *beds = bam_region2bed(in->header, region);
  } else {
    bed1_v *beds = target2bed(targets);
  }

  unsigned i;
  for (i=0; i<beds->size; ++i) {
    bed1_t *bd = ref_bed_v(beds, i);
    int nw = (bd->end - bd->beg - 1) / wsize + 1; /* number of windows */
    int copydata = calloc(nw, sizeof(int));
    
    bam_iter_t iter = bam_iter_query(idx, bd->beg+1, bd->end);
    bam1_t *b = bam_init1();
    int ret;
    while ((ret = bam_iter_read(in->x.bam, iter, b))>0) {
      bam1_core_t *c = &b->core;
      if (c->qual < conf->min_mapq) continue;
      if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->min_read_len) continue;
      if (c->flag > 0){         /* only when any flag is set */
        if (conf->filter_secondary && c->flag & BAM_FSECONDARY) continue;
        if (conf->filter_duplicate && c->flag & BAM_FDUP) continue;
        if (conf->filter_ppair && !(c->flag & BAM_FPROPER_PAIR)) continue;
        if (conf->filter_qcfail && c->flag & BAM_FQCFAIL) continue;
      }

      uint8_t *nm = bam_aux_get(b, "NM");
      if (nm && bam_aux2i(nm)>conf->max_nm) continue;
      uint32_t cnt_ret = cnt_retention(rs, b, bsstrand);
      if (cnt_ret > conf->max_retention) continue;

      uint32_t rpos = c->pos+1, qpos = 0;
      for (i=0; i<c->n_cigar; ++i) {
        uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
        uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
        switch(op) {
        case BAM_CMATCH:
          for (j=0; j<oplen; ++j)
            copydata[rpos+j / wsize]++;
          rpos += oplen;
          qpos += oplen;
          break;
        case BAM_CINS:
          qpos += oplen;
          break;
        case BAM_CDEL:
          rpos += oplen;
          break;
        case BAM_CSOFT_CLIP:
          qpos += oplen;
          break;
        case BAM_CHARD_CLIP:
          qpos += oplen;
          break;
        default:
          fprintf(stderr, "Unknown cigar, %u\n", op);
          abort();
        }
      }
    }

  
  return 0;
}
