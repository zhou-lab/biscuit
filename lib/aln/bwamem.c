#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/* Theory on probability and scoring *ungapped* alignment
 *
 * s'(a,b) = log[P(b|a)/P(b)] = log[4P(b|a)], assuming uniform base distribution
 * s'(a,a) = log(4), s'(a,b) = log(4e/3), where e is the error rate
 *
 * Scale s'(a,b) to s(a,a) s.t. s(a,a)=x. Then s(a,b) = x*s'(a,b)/log(4), or conversely: s'(a,b)=s(a,b)*log(4)/x
 *
 * If the matching score is x and mismatch penalty is -y, we can compute error rate e:
 *   e = .75 * exp[-log(4) * y/x]
 *
 * log P(seq) = \sum_i log P(b_i|a_i) = \sum_i {s'(a,b) - log(4)}
 *   = \sum_i { s(a,b)*log(4)/x - log(4) } = log(4) * (S/x - l)
 *
 * where S=\sum_i s(a,b) is the alignment score. Converting to the phred scale:
 *   Q(seq) = -10/log(10) * log P(seq) = 10*log(4)/log(10) * (l - S/x) = 6.02 * (l - S/x)
 *
 *
 * Gap open (zero gap): q' = log[P(gap-open)], r' = log[P(gap-ext)] (see Durbin et al. (1998) Section 4.1)
 * Then q = x*log[P(gap-open)]/log(4), r = x*log[P(gap-ext)]/log(4)
 *
 * When there are gaps, l should be the length of alignment matches (i.e. the M operator in CIGAR)
 */

/* static const bntseq_t *global_bns = 0; // for debugging only */

mem_opt_t *mem_opt_init() {
  mem_opt_t *o;
  o = calloc(1, sizeof(mem_opt_t));
  o->flag = 0;
  o->a = 1;
  o->b = 2;                     /* WZBS */
  /* o->b = 4; */
  o->o_del = o->o_ins = 6;
  o->e_del = o->e_ins = 1;
  o->w = 100;
  o->T = 30;
  o->zdrop = 100;
  o->pen_unpaired = 17;
  /* o->pen_clip5 = o->pen_clip3 = 5; */
  o->pen_clip5 = o->pen_clip3 = 10; /* WZBS */
  o->max_mem_intv = 20;
  o->min_seed_len = 19;
  o->split_width = 10;
  o->max_occ = 500;
  o->max_chain_gap = 10000;
  o->max_ins = 10000;
  o->mask_level = 0.50;
  o->drop_ratio = 0.50;
  o->XA_drop_ratio = 0.80;
  o->split_factor = 1.5;
  o->chunk_size = 10000000;
  o->n_threads = 1;
  o->max_XA_hits = 5;
  o->max_XA_hits_alt = 200;
  o->max_matesw = 50;
  o->mask_level_redun = 0.95;
  o->min_chain_weight = 0;
  o->max_chain_extend = 1<<30;
  o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
  o->bsstrand = 0;
  o->parent = 0;
  bwa_fill_scmat(o->a, o->b, o->mat);
  /* WZBS */
  bwa_fill_scmat_ct(o->a, o->b, o->ctmat);
  bwa_fill_scmat_ga(o->a, o->b, o->gamat);
  return o;
}

/*****************************
 * Basic hit->SAM conversion *
 *****************************/

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

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_) {
  int l_name;
  uint32_t i;
  mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert
  // "which" is 0 for primary alignment, >0 for supplementary alignment

  if (m_) mtmp = *m_, m = &mtmp;
  // set flag
  p->flag |= m? 0x1 : 0; // is paired in sequencing
  p->flag |= p->rid < 0? 0x4 : 0; // is mapped
  p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
  if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
    p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
  if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
    m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
  p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
  p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

  // print up to CIGAR
  l_name = strlen(s->name);
  ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
  kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
  kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
  if (p->rid >= 0) { // with coordinate
    kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
    kputl(p->pos + 1, str); kputc('\t', str); // POS
    kputw(p->mapq, str); kputc('\t', str); // MAPQ
    if (p->n_cigar) { // aligned
      for (i = 0; i < p->n_cigar; ++i) {
        int c = p->cigar[i]&0xf;
        if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
          c = which? 4 : 3; // use hard clipping for supplementary alignments
        kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
      }
    } else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
  } else kputsn("*\t0\t0\t*", 7, str); // without coordinte
  kputc('\t', str);

  // print the mate position if applicable
  if (m && m->rid >= 0) {
    if (p->rid == m->rid) kputc('=', str);
    else kputs(bns->anns[m->rid].name, str);
    kputc('\t', str);
    kputl(m->pos + 1, str); kputc('\t', str);
    if (p->rid == m->rid) {
      int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
      int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
      if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
      else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
    } else kputc('0', str);
  } else kputsn("*\t0\t0", 5, str);
  kputc('\t', str);

  // print SEQ and QUAL
  if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
    kputsn("*\t*", 3, str);
  } else if (!p->is_rev) { // the forward strand
    int i, qb = 0, qe = s->l_seq;
    if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
      if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
      if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
    }
    ks_resize(str, str->l + (qe - qb) + 1);
    for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
    kputc('\t', str);
    if (s->qual) { // printf qual
      ks_resize(str, str->l + (qe - qb) + 1);
      for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
      str->s[str->l] = 0;
    } else kputc('*', str);
  } else { // the reverse strand
    int i, qb = 0, qe = s->l_seq;
    if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
      if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
      if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
    }
    ks_resize(str, str->l + (qe - qb) + 1);
    for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
    kputc('\t', str);
    if (s->qual) { // printf qual
      ks_resize(str, str->l + (qe - qb) + 1);
      for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
      str->s[str->l] = 0;
    } else kputc('*', str);
  }

  // print optional tags
  if (p->n_cigar) {
    kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
    kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
    kputsn("\tZC:i:", 6, str); kputw(p->ZC, str);
    kputsn("\tZR:i:", 6, str); kputw(p->ZR, str);
  }
  if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
  if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
  if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
  if (!(p->flag & 0x100)) { // not multi-hit
    for (i = 0; i < n; ++i)
      if (i != which && !(list[i].flag&0x100)) break;
    if (i < n) { // there are other primary hits; output them
      kputsn("\tSA:Z:", 6, str);
      for (i = 0; i < n; ++i) {
        const mem_aln_t *r = &list[i];
        int k;
        if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
        kputs(bns->anns[r->rid].name, str); kputc(',', str);
        kputl(r->pos+1, str); kputc(',', str);
        kputc("+-"[r->is_rev], str); kputc(',', str);
        for (k = 0; k < r->n_cigar; ++k) {
          kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
        }
        kputc(',', str); kputw(r->mapq, str);
        kputc(',', str); kputw(r->NM, str);
        kputc(';', str);
      }
    }
    if (p->alt_sc > 0)
      ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
  }
  if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
  if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
  if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
    int tmp;
    kputsn("\tXR:Z:", 6, str);
    tmp = str->l;
    kputs(bns->anns[p->rid].anno, str);
    for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
      if (str->s[i] == '\t') str->s[i] = ' ';
  }
  /* WZBS */
  kputsn("\tYD:A:", 6, str);    /* whether it's BSW or BSC */
  if (p->bss < 0) kputc('u', str);
  else kputc("fr"[p->bss], str);

  kputc('\n', str);
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a)
{
  int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
  double identity;
  sub = a->csub > sub? a->csub : sub;
  if (sub >= a->score) return 0;
  l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
  identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
  if (a->score == 0) {
    mapq = 0;
  } else if (opt->mapQ_coef_len > 0) {
    double tmp;
    tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
    tmp *= identity * identity;
    mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
  } else {
    mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
    mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
  }
  if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
  if (mapq > 60) mapq = 60;
  if (mapq < 0) mapq = 0;
  mapq = (int)(mapq * (1. - a->frac_rep) + .499);
  return mapq;
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible

uint8_t *bseq_bsconvert(bseq1_t *s, uint8_t parent) {
  if (s->bisseq[parent]) return s->bisseq[parent];
  uint32_t i;
  if (parent) {
    s->bisseq[1] = calloc(s->l_seq, sizeof(uint8_t));
    for (i=0; i< (unsigned) s->l_seq; ++i) {
      if (s->seq[i] == 1) s->bisseq[1][i] = 3;
      else s->bisseq[1][i] = s->seq[i];
    }
    return s->bisseq[1];
  } else {
    s->bisseq[0] = calloc(s->l_seq, sizeof(uint8_t));
    for (i=0; i< (unsigned) s->l_seq; ++i) {
      if (s->seq[i] == 2) s->bisseq[0][i] = 0;
      else s->bisseq[0][i] = s->seq[i];
    }
    return s->bisseq[0];
  }
}

/**
 * @param bseq - read sequence
 * @return mem_alnreg_v* regs */
static void mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *bseq, void *buf, mem_alnreg_v *regs, uint8_t parent) {
  /* int l_seq, char *seq,  */
  int l_seq = bseq->l_seq; uint8_t *bisseq = bseq_bsconvert(bseq, parent);
  mem_chain_v chn;

  /* WZ: I think it's always 2-bit encoding */
  /* for (i = 0; i < l_seq; ++i) // convert to 2-bit encoding if we have not done so */
  /* 	seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]]; */

  /* use both bisseq and unconverted sequence here */
  chn = mem_chain(opt, bwt, bns, pac, bseq, buf, parent);
  /* filter whole chains */
  mem_chain_flt(opt, &chn);
  /* filter seeds in the chain by seed score */
  mem_flt_chained_seeds(opt, bns, pac, l_seq, bseq->seq, chn.n, chn.a, parent);
  if (bwa_verbose >= 4) mem_print_chain(bns, &chn);

  uint32_t i;
  for (i = 0; i < chn.n; ++i) {
    mem_chain_t *p = &chn.a[i];
    if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
    /* add mem_chain_t *p to mem_alnreg_v *regs */
    mem_chain2aln(opt, bns, pac, l_seq, bseq->seq, p, regs, parent);
    free(chn.a[i].seeds);
  }
  free(chn.a);
}

/**
 * make a mem_aln_t from mem_alnreg_t for sam output preparation
 * @param ar mem_alnreg_t, when 0 means unmapped
 * @return a mem_aln_t
*/
mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query_, const mem_alnreg_t *ar) {
  mem_aln_t a;
  int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1<<30), l_MD;
  int64_t pos, rb, re;
  uint8_t *query;
  uint32_t ZC, ZR;

  memset(&a, 0, sizeof(mem_aln_t));
  if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
    a.rid = -1; a.pos = -1; a.flag |= 0x4; a.bss=-1;
    return a;
  }
  qb = ar->qb, qe = ar->qe;
  rb = ar->rb, re = ar->re;
  query = malloc(l_query);
  for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
    query[i] = query_[i] < 5? query_[i] : nst_nt4_table[(int)query_[i]];
  a.mapq = ar->secondary < 0? mem_approx_mapq_se(opt, ar) : 0;
  if (ar->secondary >= 0) a.flag |= 0x100; // secondary alignment

  /* bandwidth is larger of insertion and deletion */
  tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
  w2  = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
  w2 = w2 > tmp? w2 : tmp;
  if (bwa_verbose >= 4) printf("* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
  if (w2 > opt->w) w2 = w2 < ar->w? w2 : ar->w;

  i = 0; a.cigar = 0;
  do {
    free(a.cigar);
    w2 = w2 < opt->w<<2? w2 : opt->w<<2;
    /* WZBS */
    /* ar->rid should be the real rid, see assert below */
    a.cigar = bis_bwa_gen_cigar2(ar->parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*)&query[qb], rb, re, &score, &a.n_cigar, &NM, &ZC, &ZR, ar->parent);

    if (bwa_verbose >= 4) printf("* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
    if (score == last_sc || w2 == opt->w<<2) break; // it is possible that global alignment and local alignment give different scores
    last_sc = score;
    w2 <<= 1;
  } while (++i < 3 && score < ar->truesc - opt->a);

  l_MD = strlen((char*)(a.cigar + a.n_cigar)) + 1;
  a.NM = NM;
  a.ZC = ZC;
  a.ZR = ZR;
  pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);
  a.is_rev = is_rev;
  if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
    if ((a.cigar[0]&0xf) == 2) {
      pos += a.cigar[0]>>4;
      --a.n_cigar;
      memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
    } else if ((a.cigar[a.n_cigar-1]&0xf) == 2) {
      --a.n_cigar;
      memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
    }
  }
  if (qb != 0 || qe != l_query) { // add clipping to CIGAR
    int clip5, clip3;
    clip5 = is_rev? l_query - qe : qb;
    clip3 = is_rev? qb : l_query - qe;
    a.cigar = realloc(a.cigar, 4 * (a.n_cigar + 2) + l_MD);
    if (clip5) {
      memmove(a.cigar+1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
      a.cigar[0] = clip5<<4 | 3;
      ++a.n_cigar;
    }
    if (clip3) {
      memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
      a.cigar[a.n_cigar++] = clip3<<4 | 3;
    }
  }
  a.rid = bns_pos2rid(bns, pos);
  assert(a.rid == ar->rid);
  a.pos = pos - bns->anns[a.rid].offset;
  a.score = ar->score; a.sub = ar->sub > ar->csub? ar->sub : ar->csub;
  a.is_alt = ar->is_alt; a.alt_sc = ar->alt_sc;
  a.bss = (int) ar->bss;
  free(query);
  return a;
}

static void check_paired_read_names(const char *name1, const char *name2) {
  if (strcmp(name1, name2) == 0) return;
  int l=strlen(name1);
  if (name1[l-1]=='1' && name2[l-1]=='2')
    if (strncmp(name1, name2,l-1)==0) return;
  err_fatal(__func__,"paired reads have different names: \"%s\", \"%s\"\n", name1, name2);
}

typedef struct {
  const mem_opt_t *opt;
  const bwt_t *bwt;
  const bntseq_t *bns;
  const uint8_t *pac;
  const mem_pestat_t *pes;
  bwtintv_cache_t **intv_cache;
  bseq1_t *seqs;
  mem_alnreg_v *regs;
  int64_t n_processed;
} worker_t;

/***** bisulfite adaptation *****/
/**
 * @param i i-th read is under consideration
 * @param tid thread id
 * @return w->regs[i] mem_alnreg_v*
 */
static void bis_worker1(void *data, int i, int tid)
{
  worker_t *w = (worker_t*)data;
  mem_alnreg_v *regs; const mem_opt_t *opt=w->opt;

  if (!(opt->flag&MEM_F_PE)) {	/* single-end */

    if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);

    regs = &w->regs[i]; kv_init(*regs);
    if (!(opt->parent) || !(opt->parent>>1)) /* no restriction or target daughter */
      mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i], w->intv_cache[tid], regs, 0);
    if (!(opt->parent) || opt->parent>>1) /* no restriction or target parent */
      mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i], w->intv_cache[tid], regs, 1);
    mem_merge_regions(opt, w->bns, w->pac, &w->seqs[i], regs);

  } else {			/* paired-end */

    // sanity check the read names
    check_paired_read_names(w->seqs[i<<1|0].name, w->seqs[i<<1|1].name);

    if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1 <=====\n", w->seqs[i<<1|0].name);
    regs = &w->regs[i<<1|0];
    kv_init(*regs);
    mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i<<1|0], w->intv_cache[tid], regs, 1);
    if (opt->parent)            /* align read 1 to daughter */
      mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i<<1|0], w->intv_cache[tid], regs, 0);
    mem_merge_regions(opt, w->bns, w->pac, &w->seqs[i], regs);

    if (bwa_verbose >= 4) printf("=====> Processing read '%s'/2 <=====\n", w->seqs[i<<1|1].name);
    regs = &w->regs[i<<1|1];
    kv_init(*regs);
    mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i<<1|1], w->intv_cache[tid], regs, 0);
    if (opt->parent)            /* align read 2 to parent */
      mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i<<1|1], w->intv_cache[tid], regs, 1);
    mem_merge_regions(opt, w->bns, w->pac, &w->seqs[i], regs);
  }
}

/**
 * @param i i-th read is under consideration
 * @param tid thread id
 */
static void bis_worker2(void *data, int i, int tid) {
  worker_t *w = (worker_t*)data;
  if (!(w->opt->flag&MEM_F_PE)) { /* single-end */
    if (bwa_verbose >= 4)
      printf("=====> Finalizing read '%s' <=====\n", w->seqs[i].name);

    if (w->opt->flag & MEM_F_ALN_REG) { /* output mem_alnreg_t directly */
      mem_reg2ovlp(w->opt, w->bns, &w->seqs[i], &w->regs[i]);
    } else {			/* output sam */
      mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
      mem_reg2sam(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
    }
    free(w->regs[i].a);
  } else {			/* paired-end */
    if (bwa_verbose >= 4)
      printf("=====> Finalizing read pair '%s' <=====\n", w->seqs[i<<1|0].name);

    if (!(opt->flag & MEM_F_NO_RESCUE)) 
      mem_alnreg_matesw(w->opt, w->bns, w->pac, w->pes, &w->seqs[i<<1], &w->regs[i<<1]);

    int n_pri[2];
    n_pri[0] = mem_mark_primary_se(opt, regs_pair[0].n, regs_pair[0].a, id<<1|0);
    n_pri[1] = mem_mark_primary_se(opt, regs_pair[1].n, regs_pair[1].a, id<<1|1);
    mem_sam_pe(w->opt, w->bns, w->pac, w->pes, (w->n_processed>>1) + i, &w->seqs[i<<1], &w->regs[i<<1], n_pri);
    free(w->regs[i<<1|0].a); free(w->regs[i<<1|1].a);
  }
}

/**
 * @param n: number of reads (n includes both ends for paired-end)
 * @param seqs: query sequences
 * @param pes0: paired-end statistics
 */
void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0) {

  extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
  mem_pestat_t pes[4]; int i;

  double ctime, rtime;
  ctime = cputime(); rtime = realtime();
  /* global_bns = bns;		[> get rid of this <] */

  /* initiate worker, shared across all threads */
  worker_t w;
  w.regs = malloc(n * sizeof(mem_alnreg_v));
  w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
  w.seqs = seqs; w.n_processed = n_processed;
  w.pes = &pes[0]; // isn't this shared across all threads?

  /***** Step 1: Generate mapping position *****/
  /* w.intv_cache[i] is used by thread i only */
  w.intv_cache = malloc(opt->n_threads * sizeof(bwtintv_cache_t));
  for (i = 0; i < opt->n_threads; ++i)
    w.intv_cache[i] = bwtintv_cache_init();

  kt_for(opt->n_threads, bis_worker1, &w, (opt->flag&MEM_F_PE)? n>>1 : n);

  for (i = 0; i < opt->n_threads; ++i)
    bwtintv_cache_destroy(w.intv_cache[i]);
  free(w.intv_cache);

  /***** Step 2: Obtain PE statistics *****/
  if (opt->flag & MEM_F_PE) { // infer insert sizes if not provided
    if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
    else mem_pestat(opt, bns->l_pac, n, w.regs, pes); // otherwise, infer the insert size distribution from data
  }

  /***** Step 3: Pairing and generate mapping *****/
  kt_for(opt->n_threads, bis_worker2, &w, (opt->flag&MEM_F_PE)? n>>1 : n);

  free(w.regs);

  if (bwa_verbose >= 3)
    fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}
