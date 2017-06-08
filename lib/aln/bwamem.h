#ifndef BWAMEM_H_
#define BWAMEM_H_

#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"

#define MEM_MAPQ_COEF 30.0
#define MEM_MAPQ_MAX  60

struct __smem_i;
typedef struct __smem_i smem_i;

#define MEM_F_PE        0x2
#define MEM_F_NOPAIRING 0x4
#define MEM_F_ALL       0x8
#define MEM_F_NO_MULTI  0x10
#define MEM_F_NO_RESCUE 0x20
#define MEM_F_SELF_OVLP 0x40
#define MEM_F_ALN_REG   0x80
#define MEM_F_REF_HDR	0x100
#define MEM_F_SOFTCLIP  0x200 // softclip all, by default will hardclip secondary/supplementary mapping
#define MEM_F_SMARTPE   0x400

typedef struct {
  int a, b;               // match score and mismatch penalty
  int o_del, e_del;
  int o_ins, e_ins;
  int pen_unpaired;       // phred-scaled penalty for unpaired reads
  int pen_clip5,pen_clip3;// clipping penalty. This score is not deducted from the DP score.
  int w;                  // band width
  int zdrop;              // Z-dropoff

  uint64_t max_mem_intv;

  int T;                  // output score threshold; only affecting output
  int flag;               // see MEM_F_* macros
  int min_seed_len;       // minimum seed length
  int min_chain_weight;
  uint32_t max_chain_extend;
  float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
  int split_width;        // split into a seed if its occurence is smaller than this value
  uint32_t max_occ;            // skip a seed if its occurence is larger than this value
  int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
  int n_threads;          // number of threads
  int chunk_size;         // process chunk_size-bp sequences in a batch
  float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
  float drop_ratio;       // drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
  float XA_drop_ratio;    // when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score; only effective for the XA tag
  float mask_level_redun;
  float mapQ_coef_len;
  int mapQ_coef_fac;
  int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
  int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
  int max_XA_hits, max_XA_hits_alt; // if there are max_hits or fewer, output them all
  int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
  
  /* reads can only be mapped to parent or daughter strands
   * opt->parent&1: is restricted
   * opt->parent>>1: restricted strand */
  uint8_t parent;

  /* the following would almost be impossible unless some kind of sequence context capture
   * reads can only be mapped to BSW or BSC strands
   * bsstrand&1: is restricted
   * bsstrand>>1: restricted strand */
  uint8_t bsstrand;

  /* bisulfite scoring matrix */
  int8_t ctmat[25];       /* C>T matrix */
  int8_t gamat[25];       /* G>A matrix */
} mem_opt_t;

typedef struct mem_alnreg_t mem_alnreg_t;

struct mem_alnreg_t {
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
  int pos;
  int flag;
  uint32_t is_rev:1;
  uint32_t sam_set:1;
  unsigned mapq;
  int NM;
  uint32_t ZC, ZR;
  int n_cigar;
  uint32_t *cigar;    // needs be free-ed per align
  //mem_alnreg_t *mate; // point to mate read alignment in pairing
};

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

typedef struct {
  size_t n, m;
  mem_alnreg_t *a;
  size_t n_pri; // number of regions on primary chromosomes
} mem_alnreg_v;

typedef struct {
  int low, high;   // lower and upper bounds within which a read pair is considered to be properly paired
  int failed;      // non-zero if the orientation is not supported by sufficient data
  double avg, std; // mean and stddev of the insert size distribution
} mem_pestat_t;

/* // the "finalized" version of mem_alnreg_t, it's ready for SAM output */
/* // mem_alnreg_t lacks finalized cigar, position on the chromosome and mapping quality */
/* typedef struct { // This struct is only used for the convenience of API. */
/*   int64_t pos;     // forward strand 5'-end mapping position */
/*   int rid;         // reference sequence index in bntseq_t; <0 for unmapped */
/*   int flag;        // extra flag */
/*   uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance */
/*   uint32_t ZC, ZR; */
/*   int bss;                      /\* -1: unmapped, 0: BSW, 1: BSC *\/ */
/*   int n_cigar;                  /\* can this be unsigned? number of CIGAR operations *\/ */
/*   uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234 */
/*   char *XA;        // alternative mappings */
/*   int score, sub, alt_sc; */
/* } mem_aln_t; */

typedef enum {BSS_UNSPEC, BSS_PARENT, BSS_DAUGHTER} bsstrand_t;

#ifdef __cplusplus
extern "C" {
#endif

  smem_i *smem_itr_init(const bwt_t *bwt);
  void smem_itr_destroy(smem_i *itr);
  void smem_set_query(smem_i *itr, int len, const uint8_t *query);
  void smem_config(smem_i *itr, int min_intv, int max_len, uint64_t max_intv);
  const bwtintv_v *smem_next(smem_i *itr);

  mem_opt_t *mem_opt_init(void);
  void mem_fill_scmat(int a, int b, int8_t mat[25]);

  /**
   * Align a batch of sequences and generate the alignments in the SAM format
   *
   * This routine requires $seqs[i].{l_seq,seq,name} and write $seqs[i].sam.
   * Note that $seqs[i].sam may consist of several SAM lines if the
   * corresponding sequence has multiple primary hits.
   *
   * In the paired-end mode (i.e. MEM_F_PE is set in $opt->flag), query
   * sequences must be interleaved: $n must be an even number and the 2i-th
   * sequence and the (2i+1)-th sequence constitute a read pair. In this
   * mode, there should be enough (typically >50) unique pairs for the
   * routine to infer the orientation and insert size.
   *
   * @param opt    alignment parameters
   * @param bwt    FM-index of the reference sequence
   * @param bns    Information of the reference
   * @param pac    2-bit encoded reference
   * @param n      number of query sequences
   * @param seqs   query sequences; $seqs[i].seq/sam to be modified after the call
   * @param pes0   insert-size info; if NULL, infer from data; if not NULL, it should be an array with 4 elements,
   *               corresponding to each FF, FR, RF and RR orientation. See mem_pestat() for more info.
   */
  void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0);

  static inline int get_pri_idx(double XA_drop_ratio, const mem_alnreg_t *a, int i)  {
    int k = a[i].secondary_all;
    if (k >= 0 && a[i].score >= a[k].score * XA_drop_ratio) return k;
    return -1;
  }

  /**
   * bandwidth for Smith-Waterman
   * @param a matching score
   * @param q gap opening penalty
   * @param r gap extension penalty
   **/
  static inline int infer_bw(int l1, int l2, int score, int a, int q, int r) {
    int w;
    if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
    w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
    if (w < abs(l1 - l2)) w = abs(l1 - l2);
    return w;
  }

  static inline int get_rlen(int n_cigar, const uint32_t *cigar) {
    int k, l;
    for (k = l = 0; k < n_cigar; ++k) {
      int op = cigar[k]&0xf;
      if (op == 0 || op == 2)
        l += cigar[k]>>4;
    }
    return l;
  }

  /* Only mem_align1 is kept */
  
  /* /\** */
  /*  * Find the aligned regions for one query sequence */
  /*  * */
  /*  * Note that this routine does not generate CIGAR. CIGAR should be */
  /*  * generated later by mem_reg2aln() below. */
  /*  * */
  /*  * @param opt    alignment parameters */
  /*  * @param bwt    FM-index of the reference sequence */
  /*  * @param bns    Information of the reference */
  /*  * @param pac    2-bit encoded reference */
  /*  * @param l_seq  length of query sequence */
  /*  * @param seq    query sequence */
  /*  * */
  /*  * @return       list of aligned regions. */
  /*  *\/ */
  /* mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, bsstrand_t bss); */

  /**
   * Generate CIGAR and forward-strand position from alignment region
   *
   * @param opt    alignment parameters
   * @param bns    Information of the reference
   * @param pac    2-bit encoded reference
   * @param l_seq  length of query sequence
   * @param seq    query sequence
   * @param ar     one alignment region
   *
   * @return       CIGAR, strand, mapping quality and forward-strand position
   */
  /* mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const uint8_t *seq, const mem_alnreg_t *ar); */

  /**
   * Infer the insert size distribution from interleaved alignment regions
   *
   * This function can be called after mem_align1(), as long as paired-end
   * reads are properly interleaved.
   *
   * @param opt    alignment parameters
   * @param l_pac  length of concatenated reference sequence
   * @param n      number of query sequences; must be an even number
   * @param regs   region array of size $n; 2i-th and (2i+1)-th elements constitute a pair
   * @param pes    inferred insert size distribution (output)
   */
  void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4]);

  void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, bseq1_t *s, mem_alnreg_v *a);
  /* int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]); */
  /* char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_alnreg_v *a, int l_query, const uint8_t *query); */
  void mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, const mem_alnreg_v *regs);

  int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);
  void mem_mark_primary_se(const mem_opt_t *opt, mem_alnreg_v *regs, int64_t id);
  void mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, mem_alnreg_v *regs);

  void mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], bseq1_t s[2], mem_alnreg_v regs_pair[2], int id, int *score, int *sub, int *n_sub, int z[2]);


  
  /* void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m); */
  void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *regs, const mem_alnreg_t *universal_mreg);
  void mem_reg2sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint64_t id, bseq1_t s[2], mem_alnreg_v regs_pair[2], mem_pestat_t *pes);

#ifdef __cplusplus
}
#endif

#endif
