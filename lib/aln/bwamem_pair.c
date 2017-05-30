#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "kstring.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "ksw.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif


#define MIN_RATIO     0.8
#define MIN_DIR_CNT   10
#define MIN_DIR_RATIO 0.05
#define OUTLIER_BOUND 2.0
#define MAPPING_BOUND 3.0
#define MAX_STDDEV    4.0


static int cal_sub(const mem_opt_t *opt, mem_alnreg_v *r) {
  int j;
  for (j = 1; j < r->n; ++j) { // choose unique alignment
    int b_max = r->a[j].qb > r->a[0].qb? r->a[j].qb : r->a[0].qb;
    int e_min = r->a[j].qe < r->a[0].qe? r->a[j].qe : r->a[0].qe;
    if (e_min > b_max) { // have overlap
      int min_l = r->a[j].qe - r->a[j].qb < r->a[0].qe - r->a[0].qb? r->a[j].qe - r->a[j].qb : r->a[0].qe - r->a[0].qb;
      if (e_min - b_max >= min_l * opt->mask_level) break; // significant overlap
    }
  }
  return j < r->n? r->a[j].score : opt->min_seed_len * opt->a;
}


void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m);

#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))

void check_paired_read_names(const char *name1, const char *name2) {
  if (strcmp(name1, name2) == 0) return;
  int l=strlen(name1);
  if (name1[l-1]=='1' && name2[l-1]=='2')
    if (strncmp(name1, name2,l-1)==0) return;
  err_fatal(__func__,"paired reads have different names: \"%s\", \"%s\"\n", name1, name2);
}

int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]) {

  int n = 0, i, j, z[2], o, subo, n_sub, extra_flag = 1, n_pri[2], n_aa[2];
  kstring_t str;
  mem_aln_t h[2], g[2], aa[2][2];

  str.l = str.m = 0; str.s = 0;
  memset(h, 0, sizeof(mem_aln_t) * 2);
  memset(g, 0, sizeof(mem_aln_t) * 2);
  n_aa[0] = n_aa[1] = 0;
  if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment
    mem_alnreg_v b[2];
    kv_init(b[0]); kv_init(b[1]);
    for (i = 0; i < 2; ++i)
      for (j = 0; j < a[i].n; ++j)
        if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
          kv_push(mem_alnreg_t, b[i], a[i].a[j]);
    for (i = 0; i < 2; ++i)
      for (j = 0; j < b[i].n && j < opt->max_matesw; ++j)
        n += mem_matesw(opt, bns, pac, pes, &b[i].a[j], s[!i].l_seq, (uint8_t*)s[!i].seq, &a[!i]);
    free(b[0].a); free(b[1].a);
  }
  n_pri[0] = mem_mark_primary_se(opt, a[0].n, a[0].a, id<<1|0);
  n_pri[1] = mem_mark_primary_se(opt, a[1].n, a[1].a, id<<1|1);

  if (opt->flag&MEM_F_NOPAIRING) goto no_pairing;

  /* pairing mate reads */
  if (n_pri[0] && n_pri[1] && (o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri)) > 0) {

    int is_multi[2], q_pe, score_un, q_se[2];
    char **XA[2];
    // check if an end has multiple hits even after mate-SW
    for (i = 0; i < 2; ++i) {
      for (j = 1; j < n_pri[i]; ++j)
        if (a[i].a[j].secondary < 0 && a[i].a[j].score >= opt->T) break;
      is_multi[i] = j < n_pri[i]? 1 : 0;
    }
    if (is_multi[0] || is_multi[1]) goto no_pairing; // TODO: in rare cases, the true hit may be long but with low score
    // compute mapQ for the best SE hit
    score_un = a[0].a[0].score + a[1].a[0].score - opt->pen_unpaired;
    //q_pe = o && subo < o? (int)(MEM_MAPQ_COEF * (1. - (double)subo / o) * log(a[0].a[z[0]].seedcov + a[1].a[z[1]].seedcov) + .499) : 0;
    subo = subo > score_un? subo : score_un;
    q_pe = raw_mapq(o - subo, opt->a);
    if (n_sub > 0) q_pe -= (int)(4.343 * log(n_sub+1) + .499);
    if (q_pe < 0) q_pe = 0;
    if (q_pe > 60) q_pe = 60;
    q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499);

    /* the following assumes no split hits */
    if (o > score_un) { // paired alignment is preferred
      mem_alnreg_t *c[2];
      c[0] = &a[0].a[z[0]]; c[1] = &a[1].a[z[1]];
      for (i = 0; i < 2; ++i) {
        if (c[i]->secondary >= 0)
          c[i]->sub = a[i].a[c[i]->secondary].score, c[i]->secondary = -2;
        q_se[i] = mem_approx_mapq_se(opt, c[i]);
      }
      q_se[0] = q_se[0] > q_pe? q_se[0] : q_pe < q_se[0] + 40? q_pe : q_se[0] + 40;
      q_se[1] = q_se[1] > q_pe? q_se[1] : q_pe < q_se[1] + 40? q_pe : q_se[1] + 40;
      extra_flag |= 2;
      // cap at the tandem repeat score
      q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a)? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
      q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a)? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);
    } else { // the unpaired alignment is preferred
      z[0] = z[1] = 0;
      q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
      q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
    }
    for (i = 0; i < 2; ++i) {
      int k = a[i].a[z[i]].secondary_all;
      if (k >= 0 && k < n_pri[i]) { /* switch secondary and primary if both of them are non-ALT */
        assert(a[i].a[k].secondary_all < 0);
        for (j = 0; j < a[i].n; ++j)
          if (a[i].a[j].secondary_all == k || j == k)
            a[i].a[j].secondary_all = z[i];
        a[i].a[z[i]].secondary_all = -1;
      }
    }
    if (!(opt->flag & MEM_F_ALL)) {
      for (i = 0; i < 2; ++i)
        XA[i] = mem_gen_alt(opt, bns, pac, &a[i], s[i].l_seq, s[i].seq);
    } else XA[0] = XA[1] = 0;

    /* write SAM */
    for (i = 0; i < 2; ++i) {
      h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[z[i]]);
      h[i].mapq = q_se[i];
      h[i].flag |= 0x40<<i | extra_flag;
      h[i].XA = XA[i]? XA[i][z[i]] : 0;
      aa[i][n_aa[i]++] = h[i];
      if (n_pri[i] < a[i].n) { // the read has ALT hits
        mem_alnreg_t *p = &a[i].a[n_pri[i]];
        if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
        g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p);
        g[i].flag |= 0x800 | 0x40<<i | extra_flag;
        g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
        aa[i][n_aa[i]++] = g[i];
      }
    }
    for (i = 0; i < n_aa[0]; ++i)
      mem_aln2sam(opt, bns, &str, &s[0], n_aa[0], aa[0], i, &h[1]); /* write read1 hits */
    s[0].sam = strdup(str.s); str.l = 0;
    for (i = 0; i < n_aa[1]; ++i)
      mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]); /* write read2 hits */
    s[1].sam = str.s;
    /* if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name); */
    check_paired_read_names(s[0].name, s[1].name);
    // free
    for (i = 0; i < 2; ++i) {
      free(h[i].cigar); free(g[i].cigar);
      if (XA[i] == 0) continue;
      for (j = 0; j < a[i].n; ++j) free(XA[i][j]);
      free(XA[i]);
    }
  } else goto no_pairing;
  return n;

 no_pairing:
  for (i = 0; i < 2; ++i) {
    int which = -1;
    if (a[i].n) {
      if (a[i].a[0].score >= opt->T) which = 0;
      else if (n_pri[i] < a[i].n && a[i].a[n_pri[i]].score >= opt->T)
	which = n_pri[i];
    }
    if (which >= 0) h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[which]);
    else h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, 0);
  }
  if (!(opt->flag & MEM_F_NOPAIRING) && h[0].rid == h[1].rid && h[0].rid >= 0) { // if the top hits from the two ends constitute a proper pair, flag it.
    int64_t dist;
    int d;
    d = mem_infer_dir(bns->l_pac, a[0].a[0].rb, a[1].a[0].rb, &dist);
    if (!pes[d].failed && dist >= pes[d].low && dist <= pes[d].high) extra_flag |= 2;
  }
  mem_reg2sam(opt, bns, pac, &s[0], &a[0], 0x41|extra_flag, &h[1]);
  mem_reg2sam(opt, bns, pac, &s[1], &a[1], 0x81|extra_flag, &h[0]);
  /* if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name); */
  check_paired_read_names(s[0].name, s[1].name);
  free(h[0].cigar); free(h[1].cigar);
  return n;
}
