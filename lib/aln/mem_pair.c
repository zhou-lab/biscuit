/* Pairing Paired-End Reads
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2017 Wanding.Zhou@vai.org
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

#include <math.h>
#include "mem_alnreg.h"
#include "utils.h"
#include "kvec.h"
#include "wzmisc.h"

#define MIN_RATIO     0.8
#define MIN_DIR_CNT   10
#define MIN_DIR_RATIO 0.05
#define OUTLIER_BOUND 2.0
#define MAPPING_BOUND 3.0
#define MAX_STDDEV    4.0

static int cal_sub(const mem_opt_t *opt, mem_alnreg_v *regs) {

  mem_alnreg_t *best = &regs->a[0];

  unsigned j; mem_alnreg_t *p;
  for (j = 1; j < regs->n; ++j) { // choose unique alignment
    p = &regs->a[j];
    int b_max = p->qb > best->qb ? p->qb : best->qb;
    int e_min = p->qe < best->qe? p->qe : best->qe;
    if (e_min > b_max) { // have overlap
      int min_l = p->qe - p->qb < best->qe - best->qb? p->qe - p->qb : best->qe - best->qb;
      if (e_min - b_max >= min_l * opt->mask_level) break; // significant overlap
    }
  }
  return j < regs->n ? p->score : opt->min_seed_len * opt->a;
}

typedef struct { size_t n, m; int64_t *a; } int64_v;
mem_pestat_t mem_pestat(const mem_opt_t *opt, const bntseq_t *bns, int n, const mem_alnreg_v *regs_pairs) {

  int64_v isize = {0,0,0};

  /* infer isize distribution based on the first reg from the two reads */
  int i;
  for (i = 0; i < n>>1; ++i) {
    int64_t is; // insert size
    mem_alnreg_v *regs_pair[2];
    regs_pair[0] = (mem_alnreg_v*)&regs_pairs[i<<1|0];
    regs_pair[1] = (mem_alnreg_v*)&regs_pairs[i<<1|1];

    // skip if no mapping
    if (regs_pair[0]->n == 0 || regs_pair[1]->n == 0) continue;

    mem_alnreg_t *best0 = &regs_pair[0]->a[0];
    mem_alnreg_t *best1 = &regs_pair[1]->a[0];

    // skip if sub-optimal is too close to optimal
    if (cal_sub(opt, regs_pair[0]) > MIN_RATIO * best0->score) continue;
    if (cal_sub(opt, regs_pair[1]) > MIN_RATIO * best1->score) continue;

    // skip if on different chromosome
    if (best0->rid != best1->rid) continue;
    
    // skip if on different bisulfite converted strands
    if (best0->bss != best1->bss) continue;

    if (mem_alnreg_isize(bns, best0, best1, &is))
      if (is <= opt->max_ins) kv_push(int64_t, isize, is);
  }

  if (bwa_verbose >= 3) fprintf(stderr, "[M::%s] # candidate unique pairs: %ld\n", __func__, isize.n);

  mem_pestat_t pes; memset(&pes, 0, sizeof(mem_pestat_t));
  if (isize.n < MIN_DIR_CNT) {
    fprintf(stderr, "[M:%s] There are not enough pairs for insert size inference\n", __func__);
    free(isize.a);
    pes.failed = 1;
    return pes;
  }

  // sort
  ks_introsort_64s(isize.n, isize.a);

  int p25 = isize.a[(int)(.25 * isize.n + .499)];
  int p50 = isize.a[(int)(.50 * isize.n + .499)];
  int p75 = isize.a[(int)(.75 * isize.n + .499)];

  pes.low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
  pes.high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);

  fprintf(stderr, "[M::%s] (25, 50, 75) percentile: (%d, %d, %d)\n", __func__, p25, p50, p75);
  fprintf(stderr, "[M::%s] low and high boundaries for computing mean and std.dev: (%d, %d)\n", __func__, pes.low, pes.high);

  // average
  int x;
  for (i = x = 0, pes.avg = 0; (unsigned) i < isize.n; ++i)
    if (isize.a[i] >= pes.low && isize.a[i] <= pes.high)
      pes.avg += isize.a[i], ++x;
  pes.avg /= x;

  // std
  for (i = 0, pes.std = 0; (unsigned) i < isize.n; ++i)
    if (isize.a[i] >= pes.low && isize.a[i] <= pes.high)
      pes.std += (isize.a[i] - pes.avg) * (isize.a[i] - pes.avg);
  pes.std = sqrt(pes.std / x);

  fprintf(stderr, "[M::%s] mean and std.dev: (%.2f, %.2f)\n", __func__, pes.avg, pes.std);

  // low and high
  pes.low = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
  pes.high = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
  if (pes.low > pes.avg - MAX_STDDEV * pes.std) 
    pes.low  = (int)(pes.avg - MAX_STDDEV * pes.std + .499);
  if (pes.high < pes.avg - MAX_STDDEV * pes.std) 
    pes.high = (int)(pes.avg + MAX_STDDEV * pes.std + .499);
  // if (pes.low < 1) pes.low = 1;

  fprintf(stderr, "[M::%s] low and high boundaries for proper pairs: (%d, %d)\n", __func__, pes.low, pes.high);
  free(isize.a);

  return pes;
}

// z - index of the best pair cross regs_pair[0] and regs_pair[1]
void mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const mem_pestat_t pes, mem_alnreg_v regs_pair[2], int id, int *score, int *sub, int *n_sub, int z[2]) {

  int64_t l_pac = bns->l_pac;
  pair64_v v;
  kv_init(v);

  int i; int r; // read 1 or 2
  for (r = 0; r < 2; ++r) { // loop through read number
    for (i = 0; (unsigned) i < regs_pair[r].n_pri; ++i) {
      pair64_t key;
      mem_alnreg_t *p = &regs_pair[r].a[i];
      key.x = p->rb < l_pac ? p->rb : (l_pac<<1) - 1 - p->rb; // forward position
      /* key.x = (uint64_t)e->rid<<32 | (key.x - bns->anns[e->rid].offset); */
      /* current fix, bss is the highest bit which restrict dist, not the most efficient solution TODO */
      key.x = (uint64_t)p->bss<<63 | (uint64_t)p->rid<<32 | (key.x - bns->anns[p->rid].offset);
      key.y = (uint64_t)p->score << 32 | i << 2 | (p->rb >= l_pac)<<1 | r;
      kv_push(pair64_t, v, key);
    }
  }

  // sort by location and then ascending score
  ks_introsort_128(v.n, v.a);

  pair64_v proper_pairs;  // indices of proper pairs
  // x - merged score of the whole insert + id hash
  // y - mate index in v + read index in v
  kv_init(proper_pairs);

  // O2 solution to finding proper pairs
  int k;
  for (i = 0; (unsigned) i < v.n; ++i) {
    // v is sorted ascendingly in coordinates, going backward
    for (k = i-1; k >= 0; ++k) {
      if (v.a[i].x >> 32 != v.a[k].x >> 32) break;
      if ((int64_t) (v.a[i].x & 0xffffffffU) - (int64_t) (v.a[k].x & 0xffffffffU) > max(pes.low, pes.high)) break;

      int64_t is;
      if (mem_infer_isize(v.a[k].x, v.a[i].x, (v.a[k].y>>1)&1, (v.a[i].y>>1)&1, &is) &&
          is >= pes.low && is <= pes.high) {

        /* score of the insert by merging score of the two 
         * and the insertion properness */
        double zscore = (is - pes.avg) / pes.std;
        int _score = max(0, (int)((v.a[i].y>>32) + (v.a[k].y>>32) + .721 * log(2. * erfc(fabs(zscore) * M_SQRT1_2)) * opt->a + .499)); // .721 = 1/log(4)

        pair64_t *p = kv_pushp(pair64_t, proper_pairs);
        p->y = (uint64_t) k<<32 | i;
        p->x = (uint64_t) _score<<32 | (hash_64(p->y ^ id<<8) & 0xffffffffU);
      }
    }
  }

  if (proper_pairs.n) { // found at least one proper pair

    // sort by the insert score
    ks_introsort_128(proper_pairs.n, proper_pairs.a);

    // indices of the best pair in v
    i = proper_pairs.a[proper_pairs.n-1].y >> 32;
    k = proper_pairs.a[proper_pairs.n-1].y << 32 >> 32;

    // indices of the best pair in a
    z[v.a[i].y&1] = v.a[i].y<<32>>34; // index of the best pair read 1 in reg_pairs[0]
    z[v.a[k].y&1] = v.a[k].y<<32>>34; // index of the best pair read 2 in reg_pairs[1]

    // score of the best pairing
    *score = proper_pairs.a[proper_pairs.n-1].x >> 32;

    // score of the 2nd best pairing
    *sub = proper_pairs.n > 1 ? proper_pairs.a[proper_pairs.n-2].x>>32 : 0;

    // number of other sub-optimal pairings
    int tmp = opt->a + opt->b;
    tmp = max(tmp, opt->o_del + opt->e_del);
    tmp = max(tmp, opt->o_ins + opt->e_ins);
    for (i = (long) proper_pairs.n - 2, *n_sub = 0; i >= 0; --i)
      if (*sub - (int)(proper_pairs.a[i].x>>32) <= tmp)
        ++*n_sub;

  } else *score = 0, *sub = 0, *n_sub = 0;
  free(proper_pairs.a); free(v.a);

  return;
}



