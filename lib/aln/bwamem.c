/* bwamem algorithm
 *
 * Newly added copyright in 2022
 * Copyright (c) 2022-2023 Jacob.Morrison@vai.org
 *
 * The MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
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
#include "memchain.h"
#include "mem_alnreg.h"
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
   o->max_ins = 5000;
   o->mask_level = 0.50;
   o->drop_ratio = 0.50;
   // XA_drop_ratio: secondary hit should have score at least this fraction of the primary hit
   o->XA_drop_ratio = 0.80;
   o->split_factor = 1.5;
   o->chunk_size = 10000000;
   o->n_threads = 1;
   o->max_XA_hits = 5; // max number of primary-chr secondary hits in XA
   o->max_XA_hits_alt = 5; // max number of alt-chr secondary hits in XA
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

   // clipping
   o->clip5 = 0;
   o->clip3 = 0;
   o->min_base_qual = 0;
   return o;
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a) {
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

void bseq_bsconvert(bseq1_t *s, uint8_t parent) {
  if (s->bisseq[parent]) return;

  uint32_t i;
  if (parent) {                 // C>T strand
    s->bisseq[1] = calloc(s->l_seq, sizeof(uint8_t));
    for (i=0; i< (unsigned) s->l_seq; ++i) {
      if (s->seq[i] == 1) s->bisseq[1][i] = 3;
      else s->bisseq[1][i] = s->seq[i];
    }
  } else {                      // G>A strand
    s->bisseq[0] = calloc(s->l_seq, sizeof(uint8_t));
    for (i=0; i< (unsigned) s->l_seq; ++i) {
      if (s->seq[i] == 2) s->bisseq[0][i] = 0;
      else s->bisseq[0][i] = s->seq[i];
    }
  }
}

/**
 * @param bseq - read sequence
 * @return mem_alnreg_v* regs */
static void mem_align1_core(
   const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns,
   const uint8_t *pac, bseq1_t *bseq, void *buf, mem_alnreg_v *regs,
   uint8_t parent) {

   if (bwa_verbose >= 4) 
      printf("[%s] === Seeding %s against (parent: %u)\n", __func__, bseq->name, parent);

   bseq_bsconvert(bseq, parent); // set bseq->bisseq

   /* WZ: I think it's always 2-bit encoding */
   /* for (i = 0; i < l_seq; ++i) // convert to 2-bit encoding if we have not done so */
   /* 	seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]]; */

   /* use both bisseq and unconverted sequence here */
   mem_chain_v chns = mem_chain(opt, bwt, bns, bseq, buf, parent);
   /* filter whole chains */
   mem_chain_flt(opt, &chns);
   /* filter seeds in the chain by seed score */
   /* this is not so important for short reads */
   mem_flt_chained_seeds(opt, bns, pac, bseq, &chns, parent);

   // make sure different bisulfite strand does not interfere
   mem_chain2region(opt, bns, pac, bseq, parent, &chns, regs);
   free_mem_chain_v(chns);
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
  mem_pestat_t pes;
  bwtintv_cache_t **intv_cache;
  bseq1_t *seqs;
  mem_alnreg_v *regs;
  int64_t n_processed;
} worker_t;

/*
 * The memmem() function finds the start of the first occurrence of the
 * substring 'needle' of length 'nlen' in the memory area 'haystack' of
 * length 'hlen'.
 *
 * The return value is a pointer to the beginning of the sub-string, or
 * NULL if the substring is not found.
 */
static void *wzmemmem(const void *haystack, size_t hlen, const void *needle, size_t nlen) {
  int needle_first;
  const void *p = haystack;
  size_t plen = hlen;

  if (!nlen) return NULL;

  needle_first = *(unsigned char *)needle;

  while (plen >= nlen && (p = memchr(p, needle_first, plen - nlen + 1))) {
    if (!memcmp(p, needle, nlen))
      return (void *)p;
    
    p++;
    plen = hlen - (p - haystack);
  }

  return NULL;
}

static void read_identify_adaptor(bseq1_t *seq, uint8_t *adaptor, int l_adaptor) {
  if (adaptor == NULL) seq->l_adaptor = 0;
  else {
    uint8_t *adaptor_firstbase = wzmemmem(seq->seq, seq->l_seq, adaptor, l_adaptor);
    if (adaptor_firstbase) seq->l_adaptor = seq->l_seq-(adaptor_firstbase-seq->seq);
    else {
      int i;
      for (i=l_adaptor-1; i; --i) {
        if (memcmp(seq->seq+seq->l_seq-i, adaptor, i)==0) {
          break;
        }
      }
      seq->l_adaptor = i;
    }
  }
  // seq->l_seq -= seq->l_adaptor; // shorten the read length.
}

static void clip_read_by_quality(bseq1_t *seq, int min_base_qual) {
   if (seq->qual == NULL) return;
   for (; seq->clip5 < seq->l_seq - seq->clip3; seq->clip5++) {
      if (seq->qual[seq->clip5] >= min_base_qual + 33) break;
   }
   for (; seq->l_seq - seq->clip3 >= seq->clip5; seq->clip3++) {
      if (seq->qual[seq->l_seq - seq->clip3 - 1] >= min_base_qual + 33) break;
   }
}

static void read_clipping(bseq1_t *seq, uint8_t *adaptor, int l_adaptor, const mem_opt_t *opt) {
   // clip adaptor
   read_identify_adaptor(seq, adaptor, l_adaptor);

   // clip extra base
   seq->clip5 = opt->clip5;
   seq->clip3 = opt->clip3 + seq->l_adaptor;

   // clip by base quality
   clip_read_by_quality(seq, opt->min_base_qual);

   // adjust sequence
   seq->seq0 = seq->seq; // the original sequence start
   seq->l_seq0 = seq->l_seq; // original sequence length
   seq->seq += seq->clip5;
   seq->l_seq = seq->l_seq - seq->clip3 - seq->clip5;
   if (seq->l_seq < 0) seq->l_seq = 0;
}

/***** bisulfite adaptation *****/
/**
 * @param i i-th read is under consideration
 * @param tid thread id
 * @return w->regs[i] mem_alnreg_v*
 */
static void bis_worker1(void *data, int i, int tid) {
   
   worker_t *w = (worker_t*)data;
   mem_alnreg_v *regs; const mem_opt_t *opt=w->opt;

   if (!(opt->flag&MEM_F_PE)) {	// SE

      if (bwa_verbose >= 4)
         printf("\n=====> [%s] Processing read '%s' <=====\n",
                __func__, w->seqs[i].name);

      read_clipping(&w->seqs[i], opt->adaptor1, opt->l_adaptor1, opt);
    
      regs = &w->regs[i]; kv_init(*regs); regs->n_pri = 0;
      if (!(opt->parent&1) || // no restriction
          opt->parent>>1)     // to daughter
         mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i],
                         w->intv_cache[tid], regs, 0);
    
      if (!(opt->parent&1) || // no restriction
          !(opt->parent>>1))  // to parent
         mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i],
                         w->intv_cache[tid], regs, 1);
    
      mem_merge_regions(opt, w->bns, w->pac, &w->seqs[i], regs);

   } else {			// PE

      // sanity check the read names
      check_paired_read_names(
         w->seqs[i<<1|0].name, w->seqs[i<<1|1].name);

      read_clipping(&w->seqs[i<<1|0], opt->adaptor1, opt->l_adaptor1, opt);
      read_clipping(&w->seqs[i<<1|1], opt->adaptor2, opt->l_adaptor2, opt);
    
      if (bwa_verbose >= 4)
         printf("\n=====> [%s] Processing read '%s'/1 <=====\n",
                __func__, w->seqs[i<<1|0].name);
      
      regs = &w->regs[i<<1|0];
      kv_init(*regs); regs->n_pri = 0;
      mem_align1_core(opt, w->bwt, w->bns, w->pac, &w->seqs[i<<1|0],
                      w->intv_cache[tid], regs, 1);
      
      if (!opt->parent)   /* unrestricted: align read 1 to daughter */
         mem_align1_core(opt, w->bwt, w->bns, w->pac,
                         &w->seqs[i<<1|0], w->intv_cache[tid], regs, 0);
      
      mem_merge_regions(opt, w->bns, w->pac, &w->seqs[i<<1|0], regs);

      if (bwa_verbose >= 4)
         printf("\n=====> [%s] Processing read '%s'/2 <=====\n",
                __func__, w->seqs[i<<1|1].name);
      
      regs = &w->regs[i<<1|1];
      kv_init(*regs); regs->n_pri = 0;
      
      mem_align1_core(opt, w->bwt, w->bns, w->pac,
                      &w->seqs[i<<1|1], w->intv_cache[tid], regs, 0);
      if (!opt->parent)     /* unrestricted: align read 2 to parent */
         mem_align1_core(opt, w->bwt, w->bns, w->pac,
                         &w->seqs[i<<1|1], w->intv_cache[tid], regs, 1);
      mem_merge_regions(opt, w->bns, w->pac, &w->seqs[i<<1|1], regs);
   }
}

/**
 * @param i i-th read is under consideration
 * @param tid thread id
 */
static void bis_worker2(void *data, int i, int tid) {
  (void) tid;
  worker_t *w = (worker_t*)data;

  if (!(w->opt->flag&MEM_F_PE)) { // SE
    if (bwa_verbose >= 4)
      printf("\n=====> [%s] Finalizing SE read '%s' <=====\n",
             __func__, w->seqs[i].name);

    mem_mark_primary_se(w->opt, &w->regs[i], w->n_processed + i);
    mem_alnreg_resetFLAG(&w->regs[i]);
    mem_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i]);

    mem_alnreg_freeSAM(&w->regs[i]);
    free(w->regs[i].a);
  } else {			// PE
    if (bwa_verbose >= 4)
      printf("\n=====> [%s] Finalizing PE read '%s' <=====\n",
             __func__, w->seqs[i<<1|0].name);

    if (!(w->opt->flag & MEM_F_NO_RESCUE)) 
      mem_alnreg_matesw(w->opt, w->bns, w->pac,
                        w->pes, &w->seqs[i<<1], &w->regs[i<<1]);

    if (bwa_verbose >= 4)
       printf("\n\n====== [%s] Primary-marking read 1\n", __func__);
    
    mem_mark_primary_se(w->opt, &w->regs[i<<1|0], i<<1|0);

    if (bwa_verbose >= 4)
       printf("\n\n====== [%s] Primary-marking read 2\n", __func__);
    
    mem_mark_primary_se(w->opt, &w->regs[i<<1|1], i<<1|1);

    mem_alnreg_resetFLAG(&w->regs[i<<1|0]);
    mem_alnreg_resetFLAG(&w->regs[i<<1|1]);
    mem_reg2sam_pe(
       w->opt, w->bns, w->pac, (w->n_processed>>1) + i,
       &w->seqs[i<<1], &w->regs[i<<1], w->pes);

    mem_alnreg_freeSAM(&w->regs[i<<1|0]);
    mem_alnreg_freeSAM(&w->regs[i<<1|1]);
    free(w->regs[i<<1|0].a); free(w->regs[i<<1|1].a);
  }
}

/**
 * @param n: number of reads (n includes both ends for paired-end)
 * @param seqs: query sequences
 * @param pes0: paired-end statistics
 */
void mem_process_seqs(
   const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns,
   const uint8_t *pac, int64_t n_processed, int n,
   bseq1_t *seqs, const mem_pestat_t *pes0) {

   extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
   int i;

   double ctime, rtime;
   ctime = cputime(); rtime = realtime();
   /* global_bns = bns;		[> get rid of this <] */

   /* initiate worker, shared across all threads */
   worker_t w;
   w.regs = malloc(n * sizeof(mem_alnreg_v));
   w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
   w.seqs = seqs; w.n_processed = n_processed;
   /* w.pes = pes; // isn't this shared across all threads? */

   /***** Step 1: Generate mapping position *****/
   w.intv_cache = malloc(opt->n_threads * sizeof(bwtintv_cache_t));
   for (i = 0; i < opt->n_threads; ++i)
      w.intv_cache[i] = bwtintv_cache_init(); // w.intv_cache[i] is used by thread i only

   kt_for(opt->n_threads, bis_worker1, &w, (opt->flag&MEM_F_PE)? n>>1 : n);

   for (i = 0; i < opt->n_threads; ++i)
      bwtintv_cache_destroy(w.intv_cache[i]);
   free(w.intv_cache);

   /********************************
    * Step 2: Obtain PE statistics *
    ********************************/
   if (opt->flag & MEM_F_PE) { // infer insert sizes if not provided
      if (pes0) w.pes = *pes0;
      else w.pes = mem_pestat(opt, w.bns, n, w.regs);
   }

   /***** Step 3: Pairing and generate mapping *****/
   kt_for(opt->n_threads, bis_worker2, &w, (opt->flag&MEM_F_PE)? n>>1 : n);

   free(w.regs);

   if (bwa_verbose >= 3)
      fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}

