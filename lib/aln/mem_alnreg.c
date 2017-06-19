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

#include <limits.h>
#include "mem_alnreg.h"
#include "wzmisc.h"
#include "ksort.h"
#include "ksw.h"
#include "kstring.h"
#include "kvec.h"
#include "utils.h"

/**********************
 * Merge mem_alnreg_t *
 **********************/

// sort by 1) bss; 2) ref end;
#define alnreg_slt2(a, b) ((a).bss < (b).bss || ((a).bss == (b).bss && ((a).re < (b).re)))
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

// sort by 1) bss; 2) score; 3) ref begin; 4) query begin
#define alnreg_slt(a, b) ((a).bss < (b).bss || ((a).bss == (b).bss && ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define PATCH_MAX_R_BW 0.05f
#define PATCH_MIN_SC_RATIO 0.90f

/* opt - options
   bns - reference meta
   pac - reference
   query - read sequence
   a - reg to merge, left on query
   b - reg to merge, right on query
   previously called mem_patch_reg

   return the score for concatenating mem_alnreg_t *a and mem_alnreg_t *b */
static int mem_test_reg_concatenation(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b, int *_w) {

  if (bns == 0 || pac == 0 || query == 0) return 0;

  assert(a->rid == b->rid && a->rb <= b->rb);

  if (a->rb < bns->l_pac && b->rb >= bns->l_pac) return 0; // on different strands

  if (a->qb >= b->qb || a->qe >= b->qe || a->re >= b->re) return 0; // not colinear

  // required bandwidth
  int w = (a->re - b->rb) - (a->qe - b->qb);
  w = w > 0 ? w : -w; // l = abs(l)

  // relative bandwidth
  double r = (double)(a->re - b->rb) / (b->re - a->rb) - (double)(a->qe - b->qb) / (b->qe - a->qb);
  r = r > 0.? r : -r; // r = fabs(r)

  if (bwa_verbose >= 4) printf("* potential hit merge between [%d,%d)<=>[%ld,%ld) and [%d,%d)<=>[%ld,%ld), @ %s; w=%d, r=%.4g\n", a->qb, a->qe, (long)a->rb, (long)a->re, b->qb, b->qe, (long)b->rb, (long)b->re, bns->anns[a->rid].name, w, r);

  if (a->re < b->rb || a->qe < b->qb) { // no overlap on query or on ref
    if (w > opt->w<<1 || r >= PATCH_MAX_R_BW) return 0; // the bandwidth or the relative bandwidth is too large
  } else if (w > opt->w<<2 || r >= PATCH_MAX_R_BW*2) return 0; // more permissive if overlapping on both ref and query

  // global alignment
  w += a->w + b->w;
  w = min(w, opt->w<<2);

  if (bwa_verbose >= 4) printf("* test potential hit merge with global alignment; w=%d\n", w);

  int score;
  bis_bwa_gen_cigar2(a->parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, pac, b->qe - a->qb, query + a->qb, a->rb, b->re, &score, 0, 0, 0, 0, a->parent);

  // predicted score from query
  int q_s = (int)((double)(b->qe - a->qb) / ((b->qe - b->qb) + (a->qe - a->qb)) * (b->score + a->score) + .499);
  // predicted score from ref
  int r_s = (int)((double)(b->re - a->rb) / ((b->re - b->rb) + (a->re - a->rb)) * (b->score + a->score) + .499);

  if (bwa_verbose >= 4) printf("* score=%d;(%d,%d)\n", score, q_s, r_s);

  if ((double) score / max(q_s, r_s) < PATCH_MIN_SC_RATIO) return 0;

  *_w = w;

  return score;
}

/* sort deduplicate mem_alnreg_v 
   Note when used with bns==pac==query==0, there is no merge */
void mem_sort_deduplicate(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, mem_alnreg_v *regs) {

  if (regs->n <= 1)
    return;

  // sort by the END position, not START!
  ks_introsort(mem_ars2, regs->n, regs->a);

  int i;
  for (i = 0; (unsigned) i < regs->n; ++i)
    regs->a[i].n_comp = 1;

  for (i = 1; (unsigned) i < regs->n; ++i) {

    mem_alnreg_t *p = regs->a+i;

    // compare with all previous chains with the same bss (bisulfite conversion strand),
    // rid (chromosome) and with opt->max_chain_gap distance from p
    int j;
    for (j = i - 1; j >= 0 &&
        p->bss == regs->a[j].bss && 
        p->rid == regs->a[j].rid && 
        p->rb < regs->a[j].re + opt->max_chain_gap; --j) {

      mem_alnreg_t *q = regs->a+j;

      if (q->qe == q->qb) continue; // q was excluded

      int64_t or = q->re - p->rb; // overlap length on the reference
      int64_t oq = q->qb < p->qb ? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
      int64_t mr = min(q->re - q->rb, p->re - p->rb);   // min ref len in alignment
      int64_t mq = min(q->qe - q->qb, p->qe - p->qb);   // min qry len in alignment

      int score, w;
      if (or > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
        if (p->score < q->score) {
          p->qe = p->qb;
          break;
        } else {
          q->qe = q->qb;
        }
      } else if (q->rb < p->rb && (score = mem_test_reg_concatenation(opt, bns, pac, query, q, p, &w)) > 0) {
        // merge q into p
        p->n_comp += q->n_comp + 1;
        p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
        p->sub = max(p->sub, q->sub);
        p->csub = max(p->csub, q->csub);
        p->truesc = p->score = score;
        // reset p's begin
        p->qb = q->qb;
        p->rb = q->rb;
        p->w = w;
        // mark obsolete q
        q->qb = q->qe;
      }
    }
  }

  int m;
  for (i = 0, m = 0; (unsigned) i < regs->n; ++i) { // exclude obsolete hits
    if (regs->a[i].qe > regs->a[i].qb) {
      if (m != i) regs->a[m++] = regs->a[i];
      else ++m;
    }
  }
  regs->n = m;

  // mark obsolete continguous identical hits (same score, same starting location)
  ks_introsort(mem_ars, regs->n, regs->a);
  for (i = 1; (unsigned) i < regs->n; ++i) {
    if (regs->a[i].score == regs->a[i-1].score && regs->a[i].rb == regs->a[i-1].rb && regs->a[i].qb == regs->a[i-1].qb)
      regs->a[i].qe = regs->a[i].qb;
  }

  for (i = 1, m = 1; (unsigned) i < regs->n; ++i) { // exclude identical hits
    if (regs->a[i].qe > regs->a[i].qb) {
      if (m != i) regs->a[m++] = regs->a[i];
      else ++m;
    }
  }
  regs->n = m;

  return;
}

// remove the first chain, what's the purpose?
static void mem_test_and_remove_exact(const mem_opt_t *opt, mem_alnreg_v *regs, int qlen) {
  if (!(opt->flag & MEM_F_SELF_OVLP) || regs->n == 0 || regs->a[0].truesc != qlen * opt->a) 
    return;
  memmove(regs->a, regs->a + 1, (regs->n - 1) * sizeof(mem_alnreg_t));
  regs->n--;
  return;
}

/* Merge Aligned Regions
 * previously called mem_merge_reg1 */
void mem_merge_regions(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *bseq, mem_alnreg_v *regs) {

  uint32_t i;
  mem_sort_deduplicate(opt, bns, pac, bseq->seq, regs);

  if (opt->flag & MEM_F_SELF_OVLP)
    mem_test_and_remove_exact(opt, regs, bseq->l_seq);

  if (bwa_verbose >= 4) {
    err_printf("* %ld regions remain after merging duplicated regions\n", regs->n);
    mem_print_regions(bns, regs);
  }

  /* region is on ALT chromosomes */
  for (i = 0; i < regs->n; ++i) {
    mem_alnreg_t *p = &regs->a[i];
    if (p->rid >= 0 && bns->anns[p->rid].is_alt)
      p->is_alt = 1;
  }
}



/*********************************
 * Mark primary region.
 * Two rounds of labeling were performed. 
 * The first round labels reg->secondary_all w.r.t. both primary and non-primary assembly regions.
 * The second round labels reg->secondary w.r.t. only primary assembly regions.
 *********************************/

// sort by 1) score; 2) is_alt; 3) hash;
#define alnreg_hlt(a, b)  ((a).score > (b).score || ((a).score == (b).score && ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

// sort by 1) is_alt; 2) score; 3) hash;
#define alnreg_hlt2(a, b) ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash2, mem_alnreg_t, alnreg_hlt2)
typedef kvec_t(int) int_v;

// similar to the loop in mem_chain_flt()
// n_mark - the actual regions to mark
static void mem_mark_primary_se_core(const mem_opt_t *opt, int n_mark, mem_alnreg_v *regs, int_v *z) {

  // A rough estimate of the minimum score difference between primary and secondary alignment
  int tmp = opt->a + opt->b; 
  tmp = max(opt->o_del + opt->e_del, tmp);
  tmp = max(opt->o_ins + opt->e_ins, tmp);

  z->n = 0;  // the indices of all primary alignments
  kv_push(int, *z, 0);
  int i;
  for (i = 1; i < n_mark; ++i) {
    mem_alnreg_t *a = regs->a + i;

    // check if a is a subalignment of existing primary alignment
    unsigned k; 
    for (k = 0; k < z->n; ++k) {
      mem_alnreg_t *b = regs->a + z->a[k];

      int b_max = max(a->qb, b->qb);
      int e_min = min(a->qe, b->qe);
      if (e_min > b_max) { // have overlap
        int min_l = min(a->qe - a->qb, b->qe - b->qb);
        if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
          // set a as the sub-alignment of b
          if (b->sub == 0) 
            b->sub = a->score;
          if (b->score - a->score <= tmp && (b->is_alt || !a->is_alt))
            ++b->sub_n;
          break;
        }
      }
    }

    if (k == z->n) kv_push(int, *z, i); // if a is not a subalignment, consider it as a primary alignment
    else a->secondary = z->a[k]; // otherwise, the primary of a is z->a[k] (the k-th primary alignment)
  }
}

void mem_mark_primary_se(const mem_opt_t *opt, mem_alnreg_v *regs, int64_t id) {

  if (regs->n == 0) return;

  // initiate the default of secondary labels
  int i;
  for (i = regs->n_pri = 0; (unsigned) i < regs->n; ++i) {
    mem_alnreg_t *p = regs->a + i;
    p->sub = p->alt_sc = 0;
    p->secondary = -1; // secondary to none.
    p->secondary_all = -1;
    p->hash = hash_64(id+i);
    if (!p->is_alt) ++regs->n_pri;
  }

  // high score region comes first
  ks_introsort(mem_ars_hash, regs->n, regs->a);

  /* First Round Marking */ 
  int_v z = {0,0,0};
  mem_mark_primary_se_core(opt, (int) regs->n, regs, &z);

  if (bwa_verbose >= 5) {
    for (i=0; i<regs->n; ++i) {
      mem_alnreg_t *p = regs->a + i;
      printf("in-marking (round1): %ld - %d\n", p->rb, p->secondary);
    }
  }

  /* set alt_sc - the score of primary mapping if it's on
   * an alternative chromosome */
  for (i = 0; (unsigned) i < regs->n; ++i) {
    mem_alnreg_t *p = regs->a + i;
    p->secondary_all = i; // keep the rank in the first round
    if (!p->is_alt && p->secondary >= 0 && regs->a[p->secondary].is_alt)
      p->alt_sc = regs->a[p->secondary].score;
  }

  // when there are mapping to primary chromosome and 
  // not all mappings are on primary chromosomes, remark primary mapping
  if (regs->n_pri > 0 && (unsigned) regs->n_pri < regs->n) {

    // double-dip z, expand the memory size to regs->n
    kv_resize(int, z, regs->n);

    // mapping to primary chromosomes comes first
    if (regs->n_pri > 0) 
      ks_introsort(mem_ars_hash2, regs->n, regs->a);

    // z maps rank in the 1st round to rank in the 2nd round
    for (i = 0; (unsigned) i < regs->n; ++i)
      z.a[regs->a[i].secondary_all] = i;

    // save secondary to secondary_all which is the secondary label 
    // w.r.t both primary and non-primary assemblies
    for (i = 0; (unsigned) i < regs->n; ++i) {
      if (regs->a[i].secondary >= 0) {
        regs->a[i].secondary_all = z.a[regs->a[i].secondary];
        if (regs->a[i].is_alt) regs->a[i].secondary = INT_MAX;
      } else 
        regs->a[i].secondary_all = -1;
    }

    /************************/
    /* Second Round Marking */
    /************************/
    // mark primary for hits mapped to the primary assembly only
    if (regs->n_pri > 0) {
      for (i = 0; (unsigned) i < regs->n_pri; ++i) {
        regs->a[i].sub = 0;
        regs->a[i].secondary = -1;
      }
      mem_mark_primary_se_core(opt, regs->n_pri, regs, &z);
    }
  } else {
    for (i = 0; (unsigned) i < regs->n; ++i)
      regs->a[i].secondary_all = regs->a[i].secondary;
  }


  if (bwa_verbose >= 5) {
    for (i=0; i<regs->n; ++i) {
      mem_alnreg_t *p = regs->a + i;
      printf("in-marking (round2): %ld - %d\n", p->rb, p->secondary);
    }
  }
  free(z.a);
  return;
}

/***************
 * Mate rescue *
 ***************/

/* try adding a properly-positioned mate alignment 
 * for a good-enough target alignment.
 * If success, add to mate alignments.
 * This function SW-aligns the mate sequence, and can be slow */
// regs  - target region
// l_ms  - length of mate sequence
// ms    - mate sequence (should be const?)
// mregs - mate regions
// aka mem_matesw
static void mem_alnreg_matesw_core(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes, const mem_alnreg_t *reg, int l_ms, uint8_t *ms, mem_alnreg_v *mregs) {
  
  int64_t l_pac = bns->l_pac;
  int i;
  for (i=0; (unsigned) i<mregs->n; ++i) { // if proper pair already exists
    int64_t is;
    if (mem_alnreg_isize(bns, reg, &mregs->a[i], &is) && is >= pes.low && is <= pes.high)
      return;
  }

  int read_is_rev = reg->rb >= l_pac;
  //int mate_is_rev = !read_is_rev;
  /* int is_larger = mate_is_rev; // whether the mate has larger coordinate */

  /* make the mate read sequence opposite to the direction of the primary read */
  uint8_t *mate_seq, *rev = 0;
  //if (!mate_is_rev) {
  rev = malloc(l_ms); // this is the reverse complement of ms
  for (i = 0; i < l_ms; ++i) rev[l_ms - 1 - i] = ms[i] < 4? 3 - ms[i] : 4;
  mate_seq = rev;
  //} else mate_seq = (uint8_t*) ms;

  /* determine reference boundary */
  int64_t rb = max(0, reg->rb + pes.low - l_ms);
  int64_t re = min(l_pac<<1, reg->rb + pes.high);

  /* ref is in the primary read's direction */
  uint8_t *ref = 0; int rid;
  if (rb < re) ref = bns_fetch_seq(bns, pac, &rb, (rb+re)>>1, &re, &rid);

  /* no funny things happening */
  if (reg->rid != rid || re - rb < opt->min_seed_len) { free(rev); free(ref); return; }

  // mate alignment, very slow
  /* bss !rev parent
   * 0   1    1
   * 1   1    0
   * 0   0    0
   * 1   0    1 **/
  uint8_t parent = reg->bss ^ (reg->rb < l_pac);
  int xtra = KSW_XSUBO | KSW_XSTART | (l_ms * opt->a < 250? KSW_XBYTE : 0) | (opt->min_seed_len * opt->a);
  kswr_t aln = ksw_align2(l_ms, mate_seq, re - rb, ref, 5, parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, xtra, 0);

    /* if (bwa_verbose >= 5) { */
    /*   printf("[%s] TryAdding matesw-ed region %ld-%ld. %d,%d\n", __func__, rb, re, l_ms, read_is_rev); */
    /*   mem_print_region1(bns, reg); */
    /*   printf("score %d\n", aln.score); */
    /* } */
  /* make mate mem_alnreg_t b */
  if (aln.score >= opt->min_seed_len && aln.qb >= 0) { // something goes wrong if aln.qb < 0

    mem_alnreg_t b; memset(&b, 0, sizeof(mem_alnreg_t));
    b.rid = reg->rid;
    b.is_alt = reg->is_alt;

    // mate strand is the opposite of read
    /* b.qb = read_is_rev ? l_ms - (aln.qe + 1) : aln.qb; */
    /* b.qe = read_is_rev ? l_ms - aln.qb : aln.qe + 1;  */
    /* b.rb = read_is_rev ? (l_pac<<1) - (rb + aln.te + 1) : (rb + aln.tb); */
    /* b.re = read_is_rev ? (l_pac<<1) - (rb + aln.tb) : (rb + aln.te + 1); */
    b.qb = l_ms - (aln.qe + 1);
    b.qe = l_ms - aln.qb;
    b.rb = (l_pac<<1) - (rb + aln.te + 1);
    b.re = (l_pac<<1) - (rb + aln.tb);
    b.score = aln.score;
    b.csub = aln.score2;
    b.secondary = -1;
    b.seedcov = min(b.re-b.rb, b.qe-b.qb) >> 1;
    b.bss = reg->bss;

    if (bwa_verbose >= 5) {
      printf("[%s] Adding matesw-ed region:\n", __func__);
      mem_print_region1(bns, &b);
      printf("[%s] for:\n", __func__);
      mem_print_region1(bns, reg);
    }

    // printf("*** %d, [%lld,%lld], %d:%d, (%lld,%lld), (%lld,%lld) == (%lld,%lld)\n", aln.score, rb, re, is_rev, is_larger, reg->rb, a->re, ma->a[0].rb, ma->a[0].re, b.rb, b.re);

    // insert b into ma s.t. ma remains sorted
    kv_push(mem_alnreg_t, *mregs, b); /* make room for a new element */
    for (i = 0; (unsigned) i < mregs->n - 1; ++i)    // find the insertion point 
      if (mregs->a[i].score < b.score) break;
    int tmp = i;
    for (i = mregs->n - 1; i > tmp; --i)  // move b s.t. ma remains sorted
      mregs->a[i] = mregs->a[i-1];
    mregs->a[i] = b;

    // sort deduplicate without merging
    mem_sort_deduplicate(opt, 0, 0, 0, mregs);
  }

  if (rev) free(rev);
  free(ref);
}


void mem_alnreg_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes, bseq1_t s[2], mem_alnreg_v regs_pair[2]) {

  // find good alignment regions
  mem_alnreg_v good_regs_pair[2];
  kv_init(good_regs_pair[0]); kv_init(good_regs_pair[1]);
  int i; unsigned j;
  for (i = 0; i < 2; ++i)
    for (j = 0; j < regs_pair[i].n; ++j) 
      if (regs_pair[i].a[j].score >= regs_pair[i].a[0].score - opt->pen_unpaired)
        kv_push(mem_alnreg_t, good_regs_pair[i], regs_pair[i].a[j]);

  // rescue mate alignment of good alignment if necessary
  for (i = 0; i < 2; ++i)
    for (j = 0; j < good_regs_pair[i].n && (int) j < opt->max_matesw; ++j)
      mem_alnreg_matesw_core(opt, bns, pac, pes, &good_regs_pair[i].a[j], s[!i].l_seq, (uint8_t*) s[!i].seq, &regs_pair[!i]);

  free(good_regs_pair[0].a); free(good_regs_pair[1].a);
}

/****************
 * Table Output *
 ****************/

/* instead of outputing sam, output tab-delimited table of aligned region
   currently for debugging only */
void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, bseq1_t *s, mem_alnreg_v *a) {
  uint32_t i;

  kstring_t str = {0,0,0};
  for (i = 0; i < a->n; ++i) {
    const mem_alnreg_t *p = &a->a[i];
    int is_rev, rid, qb = p->qb, qe = p->qe;
    int64_t pos, rb = p->rb, re = p->re;
    pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);
    rid = bns_pos2rid(bns, pos);
    assert(rid == p->rid);
    pos -= bns->anns[rid].offset;
    kputs(s->name, &str); kputc('\t', &str);
    kputw(s->l_seq, &str); kputc('\t', &str);
    if (is_rev) qb ^= qe, qe ^= qb, qb ^= qe; // swap
    kputw(qb, &str); kputc('\t', &str); kputw(qe, &str); kputc('\t', &str);
    kputs(bns->anns[rid].name, &str); kputc('\t', &str);
    kputw(bns->anns[rid].len, &str); kputc('\t', &str);
    kputw(pos, &str); kputc('\t', &str); kputw(pos + (re - rb), &str); kputc('\t', &str);
    ksprintf(&str, "%.3f", (double)p->truesc / opt->a / (qe - qb > re - rb? qe - qb : re - rb));
    kputc('\n', &str);
  }
  s->sam = str.s;
}
