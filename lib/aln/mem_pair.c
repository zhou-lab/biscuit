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

static inline int mem_infer_dir(int64_t l_pac, int64_t b1, int64_t b2, int64_t *dist) {
  int r1 = (b1 >= l_pac);
  int r2 = (b2 >= l_pac);

  // p2 is the coordinate of read 2 on the read 1 strand
  int64_t p2;
  if (r1 == r2) p2 = b2;
  else p2 = (l_pac<<1) - 1 - b2;

  *dist = p2 > b1 ? p2 - b1 : b1 - p2;
  return (r1 == r2 ? 0 : 1) ^ (p2 > b1 ? 0 : 3);
}

void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4]) {
  int i, d, max;
  uint64_v isize[4];
  memset(pes, 0, 4 * sizeof(mem_pestat_t));
  memset(isize, 0, sizeof(kvec_t(int)) * 4);

  /* infer isize distribution based on the first reg from the two reads */
  for (i = 0; i < n>>1; ++i) {
    int dir;
    int64_t is; // insert size
    mem_alnreg_v *r[2];
    r[0] = (mem_alnreg_v*)&regs[i<<1|0];
    r[1] = (mem_alnreg_v*)&regs[i<<1|1];
    if (r[0]->n == 0 || r[1]->n == 0) continue;
    if (cal_sub(opt, r[0]) > MIN_RATIO * r[0]->a[0].score) continue;
    if (cal_sub(opt, r[1]) > MIN_RATIO * r[1]->a[0].score) continue;
    if (r[0]->a[0].rid != r[1]->a[0].rid) continue; // not on the same chr
    if (r[0]->a[0].bss != r[1]->a[0].bss) continue; /* not on the same bisulfite strand */

    dir = mem_infer_dir(l_pac, r[0]->a[0].rb, r[1]->a[0].rb, &is);
    if (is && is <= opt->max_ins) kv_push(uint64_t, isize[dir], is);
  }

  if (bwa_verbose >= 3) fprintf(stderr, "[M::%s] # candidate unique pairs for (FF, FR, RF, RR): (%ld, %ld, %ld, %ld)\n", __func__, isize[0].n, isize[1].n, isize[2].n, isize[3].n);

  for (d = 0; d < 4; ++d) { // TODO: this block is nearly identical to the one in bwtsw2_pair.c. It would be better to merge these two.
    mem_pestat_t *r = &pes[d];
    uint64_v *q = &isize[d];
    int p25, p50, p75, x;
    if (q->n < MIN_DIR_CNT) {
      fprintf(stderr, "[M::%s] skip orientation %c%c as there are not enough pairs\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
      r->failed = 1;
      free(q->a);
      continue;
    } else fprintf(stderr, "[M::%s] analyzing insert size distribution for orientation %c%c...\n", __func__, "FR"[d>>1&1], "FR"[d&1]);

    // sort
    ks_introsort_64(q->n, q->a);

    p25 = q->a[(int)(.25 * q->n + .499)];
    p50 = q->a[(int)(.50 * q->n + .499)];
    p75 = q->a[(int)(.75 * q->n + .499)];
    r->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
    if (r->low < 1) r->low = 1;
    r->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
    fprintf(stderr, "[M::%s] (25, 50, 75) percentile: (%d, %d, %d)\n", __func__, p25, p50, p75);
    fprintf(stderr, "[M::%s] low and high boundaries for computing mean and std.dev: (%d, %d)\n", __func__, r->low, r->high);

    // average
    for (i = x = 0, r->avg = 0; i < q->n; ++i)
      if (q->a[i] >= r->low && q->a[i] <= r->high)
        r->avg += q->a[i], ++x;
    r->avg /= x;

    // std
    for (i = 0, r->std = 0; i < q->n; ++i)
      if (q->a[i] >= r->low && q->a[i] <= r->high)
        r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
    r->std = sqrt(r->std / x);

    fprintf(stderr, "[M::%s] mean and std.dev: (%.2f, %.2f)\n", __func__, r->avg, r->std);

    // low and high
    r->low  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
    r->high = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
    if (r->low  > r->avg - MAX_STDDEV * r->std) r->low  = (int)(r->avg - MAX_STDDEV * r->std + .499);
    if (r->high < r->avg - MAX_STDDEV * r->std) r->high = (int)(r->avg + MAX_STDDEV * r->std + .499);
    if (r->low < 1) r->low = 1;

    fprintf(stderr, "[M::%s] low and high boundaries for proper pairs: (%d, %d)\n", __func__, r->low, r->high);
    free(q->a);
  }

  for (d = 0, max = 0; d < 4; ++d)
    max = max > isize[d].n? max : isize[d].n;

  for (d = 0; d < 4; ++d) {
    if (pes[d].failed == 0 && isize[d].n < max * MIN_DIR_RATIO) {
      pes[d].failed = 1;
      fprintf(stderr, "[M::%s] skip orientation %c%c\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
    }
  }
}
// z - index of the best pair cross regs_pair[0] and regs_pair[1]
void mem_pair(const mem_opt_t *opt,
              const bntseq_t *bns, 
              const uint8_t *pac, 
              const mem_pestat_t pes[4], 
              bseq1_t s[2], 
              mem_alnreg_v regs_pair[2], 
              int id,       // ? read group ID? affect sorting of pairing
              int *score,   // score of the best pairing
              int *sub,     // score of 2nd best pairing
              int *n_sub,   // number of other suboptimal pairings (not including best and 2nd best)
              int z[2]) {

  int64_t l_pac = bns->l_pac;
  pair64_v v;
  kv_init(v);
  
  int i; int r; // read 1 or 2
  for (r = 0; r < 2; ++r) { // loop through read number
    for (i = 0; i < regs_pair[r].n_pri; ++i) {
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

  int last[4] = {-1,-1,-1,-1}; // keeps the last hit in the four directions
  int k;
  for (i = 0; i < v.n; ++i) {
    for (r = 0; r < 2; ++r) { /* mate direction */

      /* Allow only read and mate on different strands, i.e.,
       * pes[0].failed == pes[3].failed == 1 (00 and 11)
       * pes[1].failed == pes[2].failed == 0 (01 and 10) */
      int dir = r<<1 | (v.a[i].y>>1&1);     // mate_direction-read_direction
      if (pes[dir].failed) continue;        // invalid orientation

      /* mate direction and mate read number (1 or 2) */
      int which = r<<1 | ((v.a[i].y&1)^1);  // mate_direction-mate_read_number
      if (last[which] < 0) continue;           // no previous hits

      /* loop over mate read index
       * v is sorted with ascending score, so we are going backward */
      // TODO: this is a O(n^2) solution in the worst case; remember to check if this loop takes a lot of time (I doubt)
      for (k = lasty[which]; k >= 0; --k) {

        // mate_direction and mate_read_number has to be the same
        if ((v.a[k].y&3) != which) continue;

        int64_t dist = (int64_t)v.a[i].x - v.a[k].x;
        if (dist > pes[dir].high) break;
        if (dist < pes[dir].low)  continue;

        /* score of the insert by merging score of the two 
         * and the insertion properness */
        double zscore = (dist - pes[dir].avg) / pes[dir].std;
        int _score = max(0, (int)((v.a[i].y>>32) + (v.a[k].y>>32) + .721 * log(2. * erfc(fabs(zscore) * M_SQRT1_2)) * opt->a + .499)); // .721 = 1/log(4)

        pair64_t *p = kv_pushp(pair64_t, proper_pairs);
        p->y = (uint64_t)k<<32 | i;
        p->x = (uint64_t)_score<<32 | (hash_64(p->y ^ id<<8) & 0xffffffffU);
        //printf("[%lld,%lld]\t%d\tdist=%ld\n", v.a[k].x, v.a[i].x, q, (long)dist);
      }
    }
    lasty[v.a[i].y&3] = i;
  }

  if (proper_pairs.n) { // found at least one proper pair

    // sort by the insert score
    ks_introsort_128(proper_pairs.n, proper_pairs.a);

    /* indices of the best pair in v */
    i = proper_pairs.a[proper_pairs.n-1].y >> 32;
    k = proper_pairs.a[proper_pairs.n-1].y << 32 >> 32;

    /* indices of the best pair in a */
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



