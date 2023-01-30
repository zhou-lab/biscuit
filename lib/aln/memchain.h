/* The MIT License (MIT)
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

#include "bwamem.h"
#include "mem_alnreg.h"

/*******************
 * bwtintv_cache_t *
 *******************/

/* This struct is shared across reads processed in the same thread.
 * For performance's sake, it saves memory allocation.
 *
 * mem is the output from mem_collect_intv
 * _mem and tmpv are for internal use in mem_collect_intv
 * _mem is raw from bwt_smem1, before filtering by min_seed_len
 *
 * Previously called smem_aux_t in BWA code. */

typedef struct {
  bwtintv_v mem;
  bwtintv_v _mem;
  bwtintv_v *tmpv[2];
} bwtintv_cache_t;

static inline bwtintv_cache_t *bwtintv_cache_init() {
  bwtintv_cache_t *a;
  a = calloc(1, sizeof(bwtintv_cache_t));
  a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
  a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
  return a;
}

static inline void bwtintv_cache_destroy(bwtintv_cache_t *a) {
  free(a->tmpv[0]->a); free(a->tmpv[0]);
  free(a->tmpv[1]->a); free(a->tmpv[1]);
  free(a->mem.a); free(a->_mem.a);
  free(a);
}

/**************
 * mem_seed_t *
 **************
 * Each mem_seed_t built from a bwtintv_t.
 * To convert to normal coordinate:
 *  - pos = bns_depos(bns, rbeg, &is_rev)
 *  - then pos - bns->anns[p->rid].offset + 1 */
typedef struct {
  int64_t rbeg; // coordinate on forward-reverse reference
  int32_t qbeg, len;
  int score;
} mem_seed_t; // unaligned memory

typedef struct { size_t n, m; mem_seed_t *a;  } mem_seed_v;

/***************
 * mem_chain_t *
 ***************/

typedef struct {
   int first;           /* for internal use in mem_chain_flt, index of the first chain in overlap, -1 for not overlapping with any other seeds */
   int rid;
   uint32_t w:29;       /* weight, for sorting in mem_chain_flt; */
   uint32_t kept:2;     /* for internal book-keeping in mem_chain_flt. 0 (discard), 1 (be pointed by other seeds' first), 2 (large overlap), 3 (good) */
   uint32_t is_alt:1;
   float frac_rep;		   /* fraction of repeats */
   int64_t pos;
   mem_seed_v seeds;
   mem_seed_v seeds_extra;          // backup seeds
} mem_chain_t;

typedef struct { size_t n, m; mem_chain_t *a;  } mem_chain_v;


/********************************************
 * Cluster seeds into a chain (mem_chain_v).
 * Each chain contains one or more seeds
 ********************************************/
mem_chain_v mem_chain(
   const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns,
   bseq1_t *bseq, void *intv_cache, uint8_t parent);

// filter whole chain by chain weight and overlap with existing chains
void mem_chain_flt(const mem_opt_t *opt, mem_chain_v *chns);

// filter seeds in each chain by seed extension score
void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const bseq1_t *s, mem_chain_v *chns, uint8_t parent);

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chns);

void mem_chain2region(
   const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
   bseq1_t *bseq, uint8_t parent, mem_chain_v *chns, mem_alnreg_v *regs);

static inline void free_mem_chain_v(mem_chain_v chns) {
   unsigned i;
   for (i = 0; i < chns.n; ++i) {
      mem_chain_t *c = chns.a + i;
      free(c->seeds.a);
      free(c->seeds_extra.a);
   }
   free(chns.a);
}
