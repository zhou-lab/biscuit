/* The MIT License (MIT)
 *
 * Copyright (c) 2016-2017 Wanding.Zhou@vai.org
 * Copyright (c) 2007-2014 Attractive Chaos <attractor@live.co.uk>
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

#include <stdlib.h>
#include <math.h>
#include "memchain.h"
#include "ksort.h"
#include "kvec.h"
#include "utils.h"
#include "ksw.h"
#include "wzmisc.h"

/***********
 * Seeding *
 ***********/

/* Initial collection of SA intervals.
 * put seeds (bwtintv_t) into intv_cache (shared among reads for saving allocation)
 * @param opt mapping options
 * @param bwt, bwtc bwt index
 * @param len read length
 * @param seq read sequence, 0-3 for A-T
 * @return intv_cache (intv_cache->mem for a vector of bwtintv).  */
#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

static void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *seq, bwtintv_cache_t *intv_cache) {

  int k, x = 0, old_n;
  uint32_t i;
  int start_width = (opt->flag & MEM_F_SELF_OVLP)? 2 : 1;
  int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);

  bwtintv_v *_mem = &intv_cache->_mem;
  bwtintv_v *mem = &intv_cache->mem;
  bwtintv_v **tmpv = intv_cache->tmpv;

  // no seeds to begin with
  mem->n = 0;

  // first pass: find all SMEMs, only keep seeds with length >= min_seed_len
  while (x < len) {
    if (seq[x] < 4) {
      x = bwt_smem1(bwt, bwtc, len, seq, x, start_width, _mem, tmpv);
      for (i = 0; i < _mem->n; ++i)
        if ((uint32_t)_mem->a[i].info - (_mem->a[i].info>>32) >= (unsigned) opt->min_seed_len)
          kv_push(bwtintv_t, *mem, _mem->a[i]);
    } else ++x;
  }

  // second pass: find MEMs inside a long SMEM
  old_n = mem->n;
  for (k = 0; k < old_n; ++k) {
    bwtintv_t *p = &mem->a[k];
    int start = p->info>>32, end = (int32_t)p->info;
    if (end - start < split_len || p->x[2] > (unsigned) opt->split_width) continue;
    bwt_smem1(bwt, bwtc, len, seq, (start + end)>>1, p->x[2]+1, _mem, tmpv);
    for (i = 0; i < _mem->n; ++i)
      if ((uint32_t)_mem->a[i].info - (_mem->a[i].info>>32) >= (unsigned) opt->min_seed_len)
        kv_push(bwtintv_t, *mem, _mem->a[i]);
  }

  // third pass: LAST-like
  if (opt->max_mem_intv > 0) {
    x = 0;
    while (x < len) {
      if (seq[x] < 4) {
        if (1) {
          bwtintv_t m;
          x = bwt_seed_strategy1(bwt, bwtc, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
          if (m.x[2] > 0) kv_push(bwtintv_t, *mem, m);
        } else { // for now, we never come to this block which is slower
          x = bwt_smem1a(bwt, bwtc, len, seq, x, start_width, opt->max_mem_intv, _mem, tmpv);
          for (i = 0; i < _mem->n; ++i)
            kv_push(bwtintv_t, *mem, _mem->a[i]);
        }
      } else ++x;
    }
  }
  // sort
  ks_introsort(mem_intv, mem->n, mem->a);
}


/**************
 * mem_seed_t *
 **************/

/* filtering seed if it violates asymmetric scoring */
static int asymmetric_flt_seed(mem_seed_t *s, const uint8_t *pac, const bntseq_t *bns, bseq1_t *bseq) {
  int is_rev;
  bwtint_t pos = bns_depos(bns, s->rbeg, &is_rev);
  if (is_rev) pos -= s->len - 1;
  int64_t rb = s->rbeg;
  int64_t re = rb + s->len;
  int rid;
  uint8_t *ref = bns_fetch_seq(bns, pac, &rb, (rb+re)>>1, &re, &rid);
  int i;
  for (i=0; i<s->len; ++i) {
    /* filter seeds with T(ref)>C(read) or A(ref)>G(read) */
    if ((ref[i]==3&&bseq->seq[s->qbeg+i]==1) ||
        (ref[i]==0&&bseq->seq[s->qbeg+i]==2)) {
      free(ref);
      return 1;
    }
  }
  free(ref);
  return 0;
}

/***************
 * mem_chain_t *
 ***************/

/* chain weight is defined as:
 * min(seeds_base_cov_ref, seeds_base_cov_query) */
int mem_chain_weight(const mem_chain_t *c) {
  int64_t end;
  int j, w = 0, tmp;
  /* base coverage in query */
  for (j = 0, end = 0; j < c->n; ++j) {
    const mem_seed_t *s = &c->seeds[j];
    if (s->qbeg >= end) w += s->len;
    else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
    end = end > s->qbeg + s->len? end : s->qbeg + s->len;
  }
  tmp = w; w = 0;
  /* base coverage in reference */
  for (j = 0, end = 0; j < c->n; ++j) {
    const mem_seed_t *s = &c->seeds[j];
    if (s->rbeg >= end) w += s->len;
    else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
    end = end > s->rbeg + s->len? end : s->rbeg + s->len;
  }
  /* smaller of query coverage and reference coverage */
  w = w < tmp? w : tmp;
  return w < 1<<30? w : (1<<30)-1;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn) {
  uint32_t i; int j;
  for (i = 0; i < chn->n; ++i) {
    mem_chain_t *p = &chn->a[i];
    printf("* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
    for (j = 0; j < p->n; ++j) {
      bwtint_t pos;
      int is_rev;
      pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
      if (is_rev) pos -= p->seeds[j].len - 1;
      printf("\t%d;%d;%d,%ld(%s:%c%ld)",
             p->seeds[j].score, p->seeds[j].len, p->seeds[j].qbeg,
             (long)p->seeds[j].rbeg, bns->anns[p->rid].name,
             "+-"[is_rev], (long)(pos - bns->anns[p->rid].offset) + 1);
    }
    putchar('\n');
  }
}

/************
 * Chaining *
 ************/

/* test "mem_seed_t *s" against the last seed from the lower side
 * in "mem_chain_t *c". If the distance is close enough, merge s into c.
 * return 1 if the seed is merged into the chain
 *
 * Note that c->seeds were sorted by positions and lower than s. */
static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const mem_seed_t *s, int seed_rid) {

  const mem_seed_t *last = &c->seeds[c->n-1];

  // different chr; request a new chain
  if (seed_rid != c->rid) return 0;

  // contained seed (on both reference and query); do nothing, and report merged
  if (s->qbeg >= c->seeds[0].qbeg && s->qbeg + s->len <= last->qbeg + last->len &&
      s->rbeg >= c->seeds[0].rbeg && s->rbeg + s->len <= last->rbeg + last->len)
    return 1;

  // don't chain if on different strand
  if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) && s->rbeg >= l_pac) return 0;

  // grow the chain if
  // 1) position on ref and query are both within opt->max_chain_gap away from last seed
  // 2) the distances from last is similar (within opt->w) between on ref and on query
  int64_t qdist = s->qbeg - last->qbeg; // always non-negtive
  int64_t rdist = s->rbeg - last->rbeg;
  if (rdist >= 0 && qdist - rdist <= opt->w && rdist - qdist <= opt->w
      && qdist - last->len < opt->max_chain_gap && rdist - last->len < opt->max_chain_gap) {
    // grow the chain
    if (c->n == c->m) {
      c->m <<= 1;
      c->seeds = realloc(c->seeds, c->m * sizeof(mem_seed_t));
    }
    c->seeds[c->n++] = *s;
    return 1; // merged
  }
  return 0; // not merged, request to add a new chain
}


/*********************************************************
 * Cluster seeds (bwtintv_t) into a chain (mem_chain_v).
 * Each chain contains one or more seeds.
 *********************************************************/

#include "kbtree.h"
#define mem_getbss(parent, bns, rb) ((rb>bns->l_pac)==(parent)?1:0)
#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)
mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *bseq, void *intv_cache, uint8_t parent) {

  /* aux->mem -> kbtree_t(chn) *tree -> mem_chain_v chain */
  uint32_t i;
  int b, e, l_rep;
  int64_t l_pac = bns->l_pac;
  bwtintv_cache_t *_intv_cache;

  mem_chain_v chain;
  kv_init(chain);
  if (bseq->l_seq < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match

  /* B-tree of mem_chain_t */
  kbtree_t(chn) *tree;
  tree = kb_init(chn, KB_DEFAULT_SIZE);

  /* if cache is not given, create a temporary one */
  _intv_cache = intv_cache ? (bwtintv_cache_t*) intv_cache : bwtintv_cache_init();

  /* generate bwtintv_v (seeds) in _intv_cache->mem */
  mem_collect_intv(opt, &bwt[parent], &bwt[!parent], bseq->l_seq, bseq->bisseq[parent], _intv_cache);

  /* loop over mem and compute l_rep - number of repetitive seeds */
  for (i = 0, b = e = l_rep = 0; i < _intv_cache->mem.n; ++i) {
    bwtintv_t *intv = &_intv_cache->mem.a[i];
    if (intv->x[2] <= opt->max_occ) continue; // skip unique seeds
    int sb = (intv->info>>32), se = (uint32_t)intv->info;
    if (sb > e) l_rep += e - b, b = sb, e = se;
    else e = e > se ? e : se;
  }
  l_rep += e - b;

  /* cluster seeds into chains
   * find the closest chain from the lower side in kbtree_t(chn) *tree
   * if closest chain is nonexistent, then add the new seed as a new chain in the tree.
   * Note _intv_cache->mem is sorted by position. so this would work. */
  for (i = 0; i < _intv_cache->mem.n; ++i) {
    /* change bwtintv_t into mem_seed_t s */
    bwtintv_t *intv = &_intv_cache->mem.a[i];
    int step, slen = (uint32_t)intv->info - (intv->info>>32); /* seed length */
    uint32_t count; uint64_t k;
    // if (slen < opt->min_seed_len) continue;
    // ignore if too short or too repetitive

    /* if the interval is small (small p->x[2])
     * record every position in the interval
     * otherwise, record max_occ positions as seeds.
     * This is to avoid highly repetitive regions. lowering max_occ increase sensitivity. */
    step = intv->x[2] > opt->max_occ? intv->x[2] / opt->max_occ : 1;
    for (k = count = 0; k < intv->x[2] && count < opt->max_occ; k += step, ++count) {

      mem_chain_t tmp;    // the chain that hold the key for query
      mem_chain_t *lower; // lower bound interval in B-tree
      mem_chain_t *upper; // upper bound interval in B-tree

      /* this is the base coordinate in the forward-reverse reference */
      mem_seed_t s;
      s.rbeg = tmp.pos = bwt_sa(&bwt[parent], intv->x[0] + k);
      s.qbeg = intv->info>>32;
      s.score = s.len = slen;

      /* bridging multiple reference sequences or the forward-reverse boundary; TODO: split the seed; don't discard it!!! */
      int rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
      if (rid < 0) continue;

      /* filter seeds that do not conform to the assymetric scoring matrix */
      if (asymmetric_flt_seed(&s, pac, bns, bseq)) continue;

      /* force to a certain strand */
      if ((opt->bsstrand & 1) && mem_getbss(parent, bns, s.rbeg) != opt->bsstrand>>1) continue;

      int to_add = 0;
      if (kb_size(tree)) {
        kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
        // merge with the chain in the lower side, because bwtintv_v is sorted by position
        if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid)) to_add = 1;
      } else to_add = 1;

      /* new chain with one seed */
      if (to_add) {
        tmp.n = 1; tmp.m = 4;
        tmp.seeds = calloc(tmp.m, sizeof(mem_seed_t));
        tmp.seeds[0] = s;
        tmp.rid = rid;
        tmp.is_alt = !!bns->anns[rid].is_alt;
        kb_putp(chn, tree, &tmp);
      }
    }
  }
  if (intv_cache == 0) bwtintv_cache_destroy(_intv_cache);

  /* kbtree_t(chn) *tree to mem_chain_v *chain */
  kv_resize(mem_chain_t, chain, kb_size(tree));

  /* traverse tree and build mem_chain_v *chain */
#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
  __kb_traverse(mem_chain_t, tree, traverse_func);
#undef traverse_func

  for (i = 0; i < chain.n; ++i) chain.a[i].frac_rep = (float)l_rep / bseq->l_seq;
  if (bwa_verbose >= 4) printf("[%s] fraction of repetitive seeds: %.3f\n", __func__, (float)l_rep / bseq->l_seq);

  kb_destroy(chn, tree);
  return chain;
}

/*******************
 * Chain Filtering *
 *******************/

#define chn_beg(ch) ((ch).seeds->qbeg)
#define chn_end(ch) ((ch).seeds[(ch).n-1].qbeg + (ch).seeds[(ch).n-1].len)

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, mem_chain_t, flt_lt) /* sort by chain weight */

/* filter whole chain by chain weight and overlap with existing chains */
void mem_chain_flt(const mem_opt_t *opt, mem_chain_v *chns) {

  if (chns->n == 0) return; // no need to filter
  
  uint32_t i, k;
  kvec_t(int) to_keep = {0,0,0}; // this keeps int indices of the non-overlapping chains

  /* set and filter by chain weight */
  for (i = k = 0; i < chns->n; ++i) {
    mem_chain_t *c = chns->a+i;
    c->first = -1; c->kept = 0;
    c->w = mem_chain_weight(c);
    if (c->w < opt->min_chain_weight) free(c->seeds);
    else chns->a[k++] = *c;
  }
  chns->n = k;
  ks_introsort(mem_flt, chns->n, chns->a); /* sort by chain weight */

  /* pairwise chain comparisons, decide which chain to discard by overlap */
  chns->a[0].kept = 3;
  kv_push(int, to_keep, 0);	/* always include the "heaviest" chain */
  for (i = 1; i < chns->n; ++i) {
    int large_overlap = 0;
    for (k = 0; k < to_keep.n; ++k) { /* compare with all included chains for overlap */

      mem_chain_t *ci = &chns->a[i];
      mem_chain_t *ck = &chns->a[to_keep.a[k]];
      int b_max = chn_beg(*ck) > chn_beg(*ci) ? chn_beg(*ck) : chn_beg(*ci);
      int e_min = chn_end(*ck) < chn_end(*ci) ? chn_end(*ck) : chn_end(*ci);

      // have overlap; don't consider overlap where the kept chain is ALT while the current chain is primary
      if (e_min > b_max && (!ck->is_alt || ci->is_alt)) {
        int li = chn_end(*ci) - chn_beg(*ci);
        int lj = chn_end(*ck) - chn_beg(*ck);
        int min_l = li < lj ? li : lj;
        if (e_min - b_max >= min_l * opt->mask_level && // by default mask_level is 0.5
            min_l < opt->max_chain_gap) { // significant overlap
          large_overlap = 1;
          if (ck->first < 0) ck->first = i; // keep the first shadowed hit s.t. mapq can be more accurate
          if (ci->w < ck->w * opt->drop_ratio && ck->w - ci->w >= opt->min_seed_len<<1)
            break;
        }
      }
    }
    if (k == to_keep.n) {	/* no overlap with any included chain, then keep the chain */
      kv_push(int, to_keep, i);
      chns->a[i].kept = large_overlap ? 2 : 3;
    }
  }

  /* when a chain is pointed by another chain,
     the pointed chain is always kept */
  for (i = 0; i < to_keep.n; ++i) {
    mem_chain_t *c = chns->a + to_keep.a[i];
    if (c->first >= 0) chns->a[c->first].kept = 1;
  }
  free(to_keep.a);

  /* don't extend more than opt->max_chain_extend .kept=1/2 chains */
  for (i = k = 0; i < chns->n; ++i) {
    if (chns->a[i].kept == 0 || chns->a[i].kept == 3) continue;
    if (++k >= opt->max_chain_extend) break;
  }
  for (; i < chns->n; ++i)
    if (chns->a[i].kept < 3)
      chns->a[i].kept = 0;

  /* actual filtering of chains */
  for (i = k = 0; i < chns->n; ++i) {
    mem_chain_t *c = chns->a+i;
    if (c->kept == 0) free(c->seeds);
    else chns->a[k++] = *c;
  }

  chns->n = k;

  return;
}

/****************************
 * Testing Seeds in a chain *
 ****************************/

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200
#define MEM_HSP_COEF 1.1f
#define MEM_MINSC_COEF 5.5f
#define MEM_SEEDSW_COEF 0.05f

/* extend seeds by MEM_SHORT_EXT using Smith-Waterman, report extension score */
static int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_seed_t *s, uint8_t parent) {
  int qb, qe, rid;
  int64_t rb, re, mid, l_pac = bns->l_pac;
  uint8_t *rseq = 0;
  kswr_t x;

  if (s->len >= MEM_SHORT_LEN) return -1; // the seed is longer than the max-extend; no need to do SW
  qb = s->qbeg, qe = s->qbeg + s->len;
  rb = s->rbeg, re = s->rbeg + s->len;
  mid = (rb + re) >> 1;
  qb -= MEM_SHORT_EXT; qb = qb > 0? qb : 0;
  qe += MEM_SHORT_EXT; qe = qe < l_query? qe : l_query;
  rb -= MEM_SHORT_EXT; rb = rb > 0? rb : 0;
  re += MEM_SHORT_EXT; re = re < l_pac<<1? re : l_pac<<1;

  // the boundary case of crossing the two strands, cut off at strand boundaries
  if (rb < l_pac && l_pac < re) {
    if (mid < l_pac) re = l_pac;
    else rb = l_pac;
  }
  
  // the seed seems good enough; no need to do SW
  // also protect against fetching too long reference sequence
  if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN) return -1; 

  rseq = bns_fetch_seq(bns, pac, &rb, mid, &re, &rid);
  // the assymmetric Smith-Waterman
  x = ksw_align2(
      qe - qb, (uint8_t*) query + qb, 
      re - rb, rseq, 
      5, parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, 0);
  free(rseq);

  return x.score;
}

/* filter seeds in each chain by seed extension score */
void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const bseq1_t *s, mem_chain_v *chns, uint8_t parent) {
  
  int l_query = s->l_seq; uint8_t *query = s->seq;
  double min_l = opt->min_chain_weight ? MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
  if (min_l > MEM_SEEDSW_COEF * l_query) return; // don't run the following for short reads
  
  int min_HSP_score = (int)(opt->a * min_l + .499);
  unsigned u;
  int j, k;
  for (u = 0; u < chns->n; ++u) {
    mem_chain_t *c = chns->a+u;
    for (j = k = 0; j < c->n; ++j) {
      mem_seed_t *s = c->seeds+j;
      s->score = mem_seed_sw(opt, bns, pac, l_query, query, s, parent);
      if (s->score < 0 || s->score >= min_HSP_score) {
        s->score = s->score < 0 ? s->len * opt->a : s->score;
        c->seeds[k++] = *s;
      }
    }
    c->n = k;
  }
}

/*********************************
 * mem_chain_t >>>> mem_alnreg_t *
 *********************************/

// maximum length of a gap s.t. the alignment score is still positive
// it's then capped by 2*bandwidth (opt->w<<1)
static inline int cal_max_gap(const mem_opt_t *opt, int qlen) {
  int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
  int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
  int l = l_del > l_ins ? l_del : l_ins;
  l = l > 1 ? l : 1;
  return l < opt->w<<1? l : opt->w<<1;
}

// identify the reference sequence span against which a seed chain should be aligned against
static void mem_chain_reference_span(const mem_opt_t *opt, int l_query, int64_t l_pac, const mem_chain_t *c, int64_t rmax[]) {
  rmax[0] = l_pac<<1;
  rmax[1] = 0;

  int i;
  for (i=0; i<c->n; ++i) {
    const mem_seed_t *s = &c->seeds[i];
    // reference span for seed s: [b,e]
    int64_t b = s->rbeg - (s->qbeg + cal_max_gap(opt, s->qbeg));
    int64_t e = s->rbeg + s->len + ((l_query - s->qbeg - s->len) + cal_max_gap(opt, l_query - s->qbeg - s->len));

    rmax[0] = rmax[0] < b? rmax[0] : b;
    rmax[1] = rmax[1] > e? rmax[1] : e;
  }
  rmax[0] = rmax[0] > 0 ? rmax[0] : 0;
  rmax[1] = rmax[1] < l_pac<<1 ? rmax[1] : l_pac<<1;
  if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
    if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
    else rmax[0] = l_pac;
  }
}

/***********************
 * seed extension left *
 ***********************/

#define MAX_BAND_TRY  2

static void left_extend_seed_set_align_beg(
    const mem_opt_t *opt, // option
    const mem_seed_t *s,  // seed to extend
    const uint8_t *query, // read sequence
    uint8_t *rseq,        // reference sequence
    int64_t rmax[],       // start and end coordinate of reference sequence
    int *aw,              // bandwidth (output)
    int parent,           // parent or daughter strand?
    mem_alnreg_t *ar) {

  if (s->qbeg == 0) {
    ar->score = ar->truesc = s->len * opt->a, ar->qb = 0, ar->rb = s->rbeg;
    return;
  }

  // query sequence
  uint8_t *qs = malloc(s->qbeg);
  int i;
  for (i = 0; i < s->qbeg; ++i)
    qs[i] = query[s->qbeg - 1 - i];

  // reference sequence
  int64_t tmp = s->rbeg - rmax[0];
  uint8_t *rs = malloc(tmp);
  for (i = 0; i < tmp; ++i)
    rs[i] = rseq[tmp - 1 - i];

  int qle, tle, gtle, gscore;
  for (i = 0; i < MAX_BAND_TRY; ++i) {
    int prev = ar->score;
    *aw = opt->w << i;

    if (bwa_verbose >= 4) {
      int j;
      printf("*** [%s] Left ref:   ", __func__); 
      for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
      printf("*** [%s] Left query: ", __func__);
      for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
    }

    int max_off; // max off diagonal distance
    ar->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, *aw, opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off);

    if (bwa_verbose >= 4) { printf("*** [%s] Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", __func__, prev, ar->score, *aw, max_off); fflush(stdout); }

    if (ar->score == prev || max_off < (*aw>>1) + (*aw>>2)) break;
  }

  // check whether we prefer to reach the end of the query
  if (gscore <= 0 || gscore <= ar->score - opt->pen_clip5) { // local extension
    ar->qb = s->qbeg - qle;
    ar->rb = s->rbeg - tle;
    ar->truesc = ar->score;
  } else { // to-end extension
    ar->qb = 0;
    ar->rb = s->rbeg - gtle;
    ar->truesc = gscore;
  }
  free(qs); free(rs);
}

/************************
 * seed extension right *
 ************************/
static void right_extend_seed_set_align_end(
    const mem_opt_t *opt,  // option
    const mem_seed_t *s,   // seed to extend
    const uint8_t *query,  // read sequence
    int l_query,           // read length
    uint8_t *rseq,         // reference sequence
    int64_t rmax[],        // start and end coordinate of reference sequence
    int *aw,               // bandwidth (output)
    int parent,            // parent or daughter strand
    mem_alnreg_t *ar) {

  if (s->qbeg + s->len == l_query) {
    ar->qe = l_query, ar->re = s->rbeg + s->len;
    return;
  }

  int sc0 = ar->score;
  int qe = s->qbeg + s->len;
  int re = s->rbeg + s->len - rmax[0];
  assert(re >= 0);

  int qle, tle, gtle, gscore;
  int i;
  for (i = 0; i < MAX_BAND_TRY; ++i) {
    int prev = ar->score;
    *aw = opt->w << i;

    if (bwa_verbose >= 4) {
      int j;
      printf("*** [%s] Right ref:   ", __func__);
      for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
      printf("*** [%s] Right query: ", __func__);
      for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
    }

    int max_off;  // max off-diagonal distance
    ar->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, *aw, opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off);

    if (bwa_verbose >= 4) { printf("*** [%s] Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", __func__, prev, ar->score, *aw, max_off); fflush(stdout); }

    if (ar->score == prev || max_off < (*aw>>1) + (*aw>>2)) break;
  }

  // similar to the above
  if (gscore <= 0 || gscore <= ar->score - opt->pen_clip3) { // local extension
    ar->qe = qe + qle;
    ar->re = rmax[0] + re + tle;
    ar->truesc += ar->score - sc0;
  } else { // to-end extension
    ar->qe = l_query;
    ar->re = rmax[0] + re + gtle;
    ar->truesc += gscore - sc0;
  }
}


/* build mem_alnreg_v (arv) from mem_chain_t (c)
 * mem_alnreg_t is constructed through banded SW and seed extension 
 *
 * @param c mem_chain_t
 * @param opt parameter
 * @param bns reference meta
 * @param pac reference
 * @param l_query length of query
 * @param query query sequence, raw WITHOUT bisulfite conversion
 * @return mem_alnreg_v *av / reg, (the newly formed chain is pushed to av)
 */
void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *regs, uint8_t parent) {

  if (c->n == 0) return;

  // get reference sequence that alignment can happen: [rmax[0], rmax[1]]
  int64_t rmax[2];
  mem_chain_reference_span(opt, l_query, bns->l_pac, c, rmax);

  // retrieve the reference sequence
  int rid;
  uint8_t *rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
  assert(c->rid == rid);

  // sort seeds by score by score
  uint64_t *srt = calloc(c->n, sizeof(uint64_t));
  int i;
  for (i = 0; i < c->n; ++i)
    srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
  ks_introsort_64(c->n, srt);

  int k; unsigned u;
  for (k = c->n - 1; k >= 0; --k) { // loop from best scored seed to least
    const mem_seed_t *s = c->seeds + (uint32_t)srt[k];

    // test whether extension has been made before
    for (u = 0; u < regs->n; ++u) {
      mem_alnreg_t *reg = regs->a+u;
      int64_t rd;
      int qd, w, max_gap;

      // not fully contained
      if (s->rbeg < reg->rb || s->rbeg + s->len > reg->re || s->qbeg < reg->qb || s->qbeg + s->len > reg->qe) continue;

      // this seed may give a better alignment
      if (s->len - reg->seedlen0 > .1 * l_query) continue;

      // The seed is "around" a previous hit.
      // qd: distance ahead of the seed on query; rd: on reference
      qd = s->qbeg - reg->qb;
      rd = s->rbeg - reg->rb;
      max_gap = cal_max_gap(opt, min(qd, rd)); // the maximal gap allowed in regions ahead of the seed
      w = min(max_gap, reg->w);
      if (qd - rd < w && rd - qd < w) break;
      
      // similar to the previous four lines, but this time we look at the region behind
      qd = reg->qe - (s->qbeg + s->len);
      rd = reg->re - (s->rbeg + s->len);
      max_gap = cal_max_gap(opt, min(qd, rd));
      w = min(max_gap, reg->w);
      if (qd - rd < w && rd - qd < w) break;
    }

    // the seed is (almost) contained in an existing alignment; 
    // further testing is needed to confirm it is not leading to a different aln
    if (u < regs->n) {

      if (bwa_verbose >= 4)
        printf("** [%s] Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n", __func__, k, (long)s->len, (long)s->qbeg, (long)s->rbeg, regs->a[u].qb, regs->a[u].qe, (long)regs->a[u].rb, (long)regs->a[u].re);

      for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
        const mem_seed_t *t;
        if (srt[i] == 0) continue;
        t = &c->seeds[(uint32_t)srt[i]];
        if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
        if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
        if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
      }

      if (i == c->n) { // no overlapping seeds; then skip extension
        srt[k] = 0; // mark that seed extension has not been performed
        continue;
      }

      if (bwa_verbose >= 4)
        printf("** [%s] Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", __func__, k);
    }

    /**** seed extension ****/
    int aw[2]; // aw: actual bandwidth used in extension aw[0] - left extension; a[1] - right extension
    mem_alnreg_t *reg = kv_pushp(mem_alnreg_t, *regs);
    memset(reg, 0, sizeof(mem_alnreg_t)); // reset every new mem_alnreg_t
    reg->w = aw[0] = aw[1] = opt->w;
    reg->score = reg->truesc = -1;
    reg->rid = c->rid;

    if (bwa_verbose >= 4) err_printf("** ---> [%s] Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", __func__, k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);

    left_extend_seed_set_align_beg(opt, s, query, rseq, rmax, &aw[0], parent, reg);
    right_extend_seed_set_align_end(opt, s, query, l_query, rseq, rmax, &aw[1], parent, reg);

    reg->bss = mem_getbss(parent, bns, reg->rb); /* set bisulfite strand */
    reg->parent = parent;

    if (bwa_verbose >= 4) printf("*** [%s] Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", __func__, reg->qb, reg->qe, (long) reg->rb, (long)reg->re, reg->score, aw[0], aw[1]);
 
    /* compute seedcov */
    for (i = 0, reg->seedcov = 0; i < c->n; ++i) {
      const mem_seed_t *t = &c->seeds[i];
      if (t->qbeg >= reg->qb && t->qbeg + t->len <= reg->qe && t->rbeg >= reg->rb && t->rbeg + t->len <= reg->re) // seed fully contained
        reg->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
    }
    reg->w = aw[0] > aw[1]? aw[0] : aw[1];
    reg->seedlen0 = s->len; // length of best scored seed

    reg->frac_rep = c->frac_rep;
  }
  free(srt); free(rseq);
}



