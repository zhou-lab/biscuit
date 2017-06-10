#ifndef _MEM_ALNREG_H
#define _MEM_ALNREG_H

#include "bwamem.h"

/***********************
 * mem_alnreg_t
 ***********************/
// This struct holds the alignment regions, each gets turned into a SAM record.
typedef struct {
  int64_t rb, re; // [rb,re): reference sequence in the alignment (forward-reverse coordinates)
  int qb, qe;     // [qb,qe): query sequence in the alignment
  int rid;        // reference seq ID
  int score;      // best local SW score
  int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
  int sub;        // 2nd best SW score
  int alt_sc;     /* score of primary mapping in secondary mapping if that primary is on alternative chromosome, see mem_mark_primary_se */
  int csub;       // SW score of a tandem hit
  int sub_n;      // approximate number of suboptimal hits
  int w;          // actual band width used in extension
  int seedcov;    // length of regions coverged by seeds
  int secondary;  // index of the parent hit shadowing the current hit; <0 if primary, this was done only within primary assemblies
  int secondary_all; // index of the parent hit shadowing the current hit; this was done with both primary and non-primary assemblies
  int seedlen0;   // length of the starting/best-scored seed
  int n_comp:30;  // number of sub-alignments chained together
  int is_alt:2;   // reference is an alternative chromosome
  float frac_rep;
  uint64_t hash;
  uint8_t bss:1;
  uint8_t parent:1;
  uint8_t read_in_pair:1;

  // SAM meta-information, e.g., CIGAR, mapq, sub = max(sub, csub)
  int pos, flag, NM, n_cigar;
  uint32_t is_rev:1;
  uint32_t sam_set:1;
  unsigned mapq;
  uint32_t ZC, ZR;
  uint32_t *cigar;    // needs be free-ed per align
  //mem_alnreg_t *mate; // point to mate read alignment in pairing
} mem_alnreg_t;

typedef struct {
  size_t n, m;
  mem_alnreg_t *a;
  size_t n_pri; // number of regions on primary chromosomes
} mem_alnreg_v;

// 1 for proper pairing, 0 for improper pairing
static inline void mem_alnreg_infer_isize(int64_t l_pac, mem_alnreg_t *p, mem_alnreg_t *q, int *proper, int *isize) {
  int str_p = p->rb >= l_pac;
  int str_q = q->rb >= l_pac;
  if (str_p && !str_q) {
    *isize = (l_pac<<1) - 1 - p->rb - q->rb;
    *proper = 1;
    return;
  } else if (str_q && !str_p) {
    *isize = (l_pac<<1) - 1 - q->rb - p->rb;
    *proper = 1;
    return;
  } else {
    int qq = (l_pac<<1) - 1 - q->rb;
    *isize = p->rb > qq ? p->rb - qq : qq - p->rb;
    *proper = 0;
    return;
  }
}

// Merge aligned regions, aka mem_merge_reg1
void mem_merge_regions(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *bseq, mem_alnreg_v *regs);

static inline int get_pri_idx(double XA_drop_ratio, const mem_alnreg_t *a, int i)  {
  int k = a[i].secondary_all;
  if (k >= 0 && a[i].score >= a[k].score * XA_drop_ratio) return k;
  return -1;
}


mem_pestat_t mem_pestat(const mem_opt_t *opt, int n, const mem_alnreg_v *regs_pairs);

void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, bseq1_t *s, mem_alnreg_v *a);
/* int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]); */
/* char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_alnreg_v *a, int l_query, const uint8_t *query); */
void mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, const mem_alnreg_v *regs);

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);
void mem_mark_primary_se(const mem_opt_t *opt, mem_alnreg_v *regs, int64_t id);
void mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, mem_alnreg_v *regs);


/*****************************************************
 * id - read group ID? affect sorting of pairing
 *
 * output
 * score  - score of the best pairing
 * sub    - score of the 2nd best pairing
 * n_sub  - number of other sub-optimal pairings
 *          (not including best and 2nd best)
 * z[2]   - index of the best pair cross regs_pair[0] 
 *        - and regs_pair[1]
 *****************************************************/
void mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const mem_pestat_t pes, mem_alnreg_v regs_pair[2], int id, int *score, int *sub, int *n_sub, int z[2]);

/* void mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes, const mem_alnreg_t *reg, int l_ms, const uint8_t *ms, mem_alnreg_v *mregs); */
void mem_alnreg_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes, bseq1_t s[2], mem_alnreg_v regs_pair[2]);
  
/* void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m); */
void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *regs, const mem_alnreg_t *universal_mreg);
void mem_reg2sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint64_t id, bseq1_t s[2], mem_alnreg_v regs_pair[2], mem_pestat_t pes);

#endif /* _MEM_ALNREG_H */
