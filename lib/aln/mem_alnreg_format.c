/* SAM formating for mem_regaln_t
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
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
#include <stdlib.h>
#include <limits.h>
#include "mem_alnreg.h"
#include "kvec.h"
#include "kstring.h"
#include "wzmisc.h"
/* #include "sam.h" */


// set CIGAR, pos, is_rev
// mapQ and most flags (except reverse strand) are set outside this function
static void mem_alnreg_setSAM(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_t *reg) {

  if (reg->n_cigar > 0) return;
  
  // nt4 encoding
  uint8_t *query = malloc(s->l_seq); int i;
  for (i = 0; i < s->l_seq; ++i)
    query[i] = s->seq[i] < 5 ? s->seq[i] : nst_nt4_table[(int)s->seq[i]];

  // cigar and pos
  // initial bandwidth is the max of insertion and deletion calculation
  int _w1 = infer_bw(reg->qe-reg->qb, reg->re-reg->rb, reg->truesc, opt->a, opt->o_del, opt->e_del);
  int _w2 = infer_bw(reg->qe-reg->qb, reg->re-reg->rb, reg->truesc, opt->a, opt->o_ins, opt->e_ins);
  int w = max(_w1, _w2);
  if (w > opt->w) w = min(w, reg->w);

  if (bwa_verbose >= 4) {
    printf("[%s] Generate cigar for\n", __func__);
    mem_print_region1(bns, reg);
    putchar('\n');
  }

  // incrementally double bandwidth
  uint32_t *cigar = 0; int n_cigar;
  int score; int last_sc = -(1<<30);
  for (i=0; i<3; ++i, w<<=1, last_sc=score) {
    free(cigar);
    w = min(w, opt->w<<2);

    // regenerate cigar related info with new bandwidth
    cigar = bis_bwa_gen_cigar2(reg->parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, pac, reg->qe - reg->qb, (uint8_t*) &query[reg->qb], reg->rb, reg->re, &score, &n_cigar, &reg->NM, &reg->ZC, &reg->ZR, &reg->bss_u, reg->parent);

    if (bwa_verbose >= 4) printf("[%s] w=%d, global_sc=%d, local_sc=%d\n", __func__, w, score, reg->truesc);

    if (score == last_sc) break;
    if (w == opt->w << 2) break;
    if (score >= reg->truesc - opt->a) break;
  }
  int l_MD = strlen((char*) (cigar+n_cigar))+1;

  // pos and is_rev
  int _is_rev;
  int64_t rpos = bns_depos(bns, reg->rb < bns->l_pac ? reg->rb : reg->re-1, &_is_rev);
  reg->is_rev = _is_rev;
  reg->flag |= reg->is_rev ? 0x10 : 0;

  // squeeze out leading and trailing deletions
  if (n_cigar > 0) { // squeeze out leading or trailing deletions
    if ((cigar[0]&0xf) == 2) {
      rpos += cigar[0]>>4;
      --n_cigar;
      memmove(cigar, cigar + 1, n_cigar * 4 + l_MD);
    } else if ((cigar[n_cigar-1]&0xf) == 2) {
      --n_cigar;
      memmove(cigar + n_cigar, cigar + n_cigar + 1, l_MD); // MD needs to be moved accordingly
    }
  }

  // add clipping to CIGAR
  if (reg->qb != 0 || reg->qe != s->l_seq || s->clip5 || s->clip3) {
    int clip5, clip3;
    clip5 = reg->is_rev ? s->l_seq - reg->qe + s->clip3 : reg->qb + s->clip5;
    clip3 = reg->is_rev ? reg->qb + s->clip5 : s->l_seq - reg->qe + s->clip3;
    cigar = realloc(cigar, 4 * (n_cigar + 2) + l_MD);
    if (clip5) {
      memmove(cigar+1, cigar, n_cigar * 4 + l_MD); // make room for 5'-end clipping
      cigar[0] = clip5<<4 | 3;
      ++n_cigar;
    }
    if (clip3) {
      memmove(cigar + n_cigar + 1, cigar + n_cigar, l_MD); // make room for 3'-end clipping
      cigar[n_cigar++] = clip3<<4 | 3;
    }
  }
  free(query);
  reg->n_cigar = n_cigar;
  if (reg->n_cigar > 0) reg->cigar = cigar;
  else free(cigar);

  assert(bns_pos2rid(bns, rpos) == reg->rid);
  reg->pos = rpos - bns->anns[reg->rid].offset;

  return;
}

/* Generate XA, XB put to str. */
static void mem_alnreg_tagXAXB(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, const mem_alnreg_t *p0, const mem_alnreg_v *regs0, kstring_t *sam_str) {

  // no need to set XA if all alignments are output as records
  if (!regs0 || (opt->flag & MEM_F_ALL)) return;

  int cnt_pri = 0, cnt_alt= 0; // number of primary-chr, and alt-chr secondaries
  unsigned i;
  for (i=0; i<regs0->n; ++i) {
    int r = get_pri_idx(opt->XA_drop_ratio, regs0->a, i);
    if (r >= 0 && regs0->a + r == p0) { // i is the secondary to r/p0
      if (regs0->a[i].is_alt) ++cnt_alt;
      else ++cnt_pri;
    }
  }
  
  // if the number of primary-chr secondaries is higher than opt->max_XA_hits
  // or if the number of alt-chr secondaries is higher than opt->max_XA_hits_alt
  // suppress reporting secondary mapping and only output number of secondaries
  if (cnt_pri <= opt->max_XA_hits && cnt_alt <= opt->max_XA_hits_alt) {

    kstring_t str = {0,0,0};
    int n;
    for (i=0, n=0; i<regs0->n; ++i) {

      mem_alnreg_t *q = regs0->a + i;
      int r = get_pri_idx(opt->XA_drop_ratio, regs0->a, i);
      if (r < 0 || regs0->a + r != p0) continue;

      // try set cigar if haven't yet
      if (q->n_cigar == 0) {
        mem_alnreg_setSAM(opt, bns, pac, s, q);
        if (q->n_cigar == 0) continue;
      }
    
      if (n) kputc(';', &str);
      kputs(bns->anns[q->rid].name, &str);
      kputc(',', &str); 
      kputc("+-"[q->is_rev], &str);
      kputl(q->pos + 1, &str);
      kputc(',', &str);

      int k;
      for (k = 0; k < q->n_cigar; ++k) {
        kputw(q->cigar[k]>>4, &str);
        kputc("MIDSHN"[q->cigar[k]&0xf], &str);
      }
      kputc(',', &str);
      kputw(q->NM, &str);
      ++n;
    }

    if (str.l) {
      kputsn("\tXA:Z:", 6, sam_str); 
      kputs(str.s, sam_str);
    }
    free(str.s);
  }

  // XB: number of alternative alignments
  if (cnt_pri > 0 || cnt_alt > 0) { // only when there is alternative(s)
    kputsn("\tXB:Z:", 6, sam_str);
    kputw(cnt_pri, sam_str);
    kputc(',', sam_str);
    kputw(cnt_alt, sam_str);
  }
}

/* Generate SA-tag, put to str */
static void mem_alnreg_tagSA(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, const mem_alnreg_t *p0, const mem_alnreg_v *regs0, kstring_t *sam_str) {

  if (!regs0 || p0->flag & 0x100) return;

  kstring_t str = {0,0,0};
  unsigned i;
  for (i=0; i < regs0->n; ++i) {
    mem_alnreg_t *q = regs0->a + i;
    if (q == p0 || q->n_cigar == 0 || q->flag & 0x100) continue;

    // try set cigar if haven't yet
    if (q->n_cigar == 0) {
      mem_alnreg_setSAM(opt, bns, pac, s, q);
      if (q->n_cigar == 0) continue;
    }

    kputs(bns->anns[q->rid].name, &str); kputc(',', &str);
    kputl(q->pos+1, &str); kputc(',', &str);
    kputc("+-"[q->is_rev], &str); kputc(',', &str);
    int k;
    for (k = 0; k < q->n_cigar; ++k) {
      kputw(q->cigar[k]>>4, &str);
      kputc("MIDSH"[q->cigar[k]&0xf], &str);
    }
    kputc(',', &str); kputw(q->mapq, &str);
    kputc(',', &str); kputw(q->NM, &str);
    kputc(';', &str);
  }

  if (str.l) {
    kputsn("\tSA:Z:", 6, sam_str);
    kputs(str.s, sam_str);
  }
  free(str.s);
}

/*************************
 * format SAM
 *************************/
// mate is set at final stage because the mate might be asymmetric with alternative mappings
// It doesn't not change the mem_alnreg_t inputs, since one alignment may be paired 
// with multiple other alignments 
void mem_alnreg_formatSAM(
   const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
   kstring_t *str, bseq1_t *s, const mem_alnreg_t *p0, const mem_alnreg_t *m0,
   const mem_alnreg_v *regs0, int is_primary, mem_pestat_t *pes) {

   // make copies
   mem_alnreg_t p = *p0;
   mem_alnreg_t m; memset(&m, 0, sizeof(mem_alnreg_t));
   if (m0) m = *m0;

   // set mate-related flags
   p.flag |= m0 ? 0x1 : 0; // is paired in sequencing
   p.flag |= m0 && m.rid < 0 ? 0x8 : 0; // is mate mapped

   // if the mate has certain bss so should the target read
   if (m0 && m0->bss_u == 0) p.bss_u = 0;

   if (p.rid >= 0 && m0 && m.rid >= 0 && pes && is_proper_pair(bns, &p, &m, *pes)) {
      p.flag |= 2; m.flag |= 2;
   }
   // copy mate coordinate to alignment
   if (p.rid < 0 && m0 && m.rid >= 0) {
      p.rid = m.rid; 
      p.pos = m.pos;
      p.is_rev = m.is_rev; // not 100% sure if we should copy is_rev
      p.n_cigar = 0;
   }
   // copy alignment coordinates to mate
   if (m0 && (m.rid < 0) && p.rid >= 0) {  // copy alignment to mate
      m.rid = p.rid;
      m.pos = p.pos;
      m.is_rev = p.is_rev;
      m.n_cigar = 0;
   }
   p.flag |= m0 && m.is_rev ? 0x20 : 0; // is mate on the reverse strand

   // print up to CIGAR
   int l_name = strlen(s->name);
   ks_resize(str, str->l + s->l_seq0 + l_name + (s->qual ? s->l_seq0 : 0) + 20);
   kputsn(s->name, l_name, str); kputc('\t', str); // read name, QNAME
   kputw((p.flag & 0xffff) | (p.flag & 0x10000 ? 0x100 : 0), str); kputc('\t', str); // FLAG
   if (p.rid >= 0) { // with coordinate
      kputs(bns->anns[p.rid].name, str); kputc('\t', str); // reference/chromosome name, RNAME
      kputl(p.pos + 1, str); kputc('\t', str); // POS
      kputw(p.mapq, str); kputc('\t', str); // MAPQ
      if (p.n_cigar) { // CIGAR
         int i;
         for (i = 0; i < p.n_cigar; ++i) {
            int c = p.cigar[i] & 0xf;
            if (!(opt->flag & MEM_F_SOFTCLIP) && !p.is_alt && (c == 3 || c == 4))
               c = is_primary ? 3 : 4; // use hard clipping for supplementary alignments
            kputw(p.cigar[i]>>4, str); kputc("MIDSH"[c], str);
         }
      } else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
   } else kputsn("*\t0\t0\t*", 7, str); // without coordinte
   kputc('\t', str);

   // print the mate position if applicable
   if (m0 && m.rid >= 0) {
      if (p.rid == m.rid) kputc('=', str);
      else kputs(bns->anns[m.rid].name, str);
      kputc('\t', str);
      kputl(m.pos + 1, str); kputc('\t', str);
      if (p.rid == m.rid) {

         // the following calculation of insert size is different from BWA
         int64_t p0 = -1, p1 = -1;
         if (p.is_rev) p1 = p.pos + get_rlen(p.n_cigar, p.cigar) - 1;
         else p0 = p.pos;
         if (m.is_rev) p1 = m.pos + get_rlen(m.n_cigar, m.cigar) - 1;
         else p0 = m.pos;
         if (p.n_cigar > 0 && m.n_cigar > 0 && p0 >= 0 && p1 >= 0) kputl(p1-p0+1, str);
         else kputc('0', str);

         // the BWA way
         // int64_t p0 = p.pos + (p.is_rev? get_rlen(p.n_cigar, p.cigar) - 1 : 0);
         // int64_t p1 = m.pos + (m.is_rev? get_rlen(m.n_cigar, m.cigar) - 1 : 0);
         // if (m.n_cigar == 0 || p.n_cigar == 0) kputc('0', str);
         // else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
      } else kputc('0', str);
   } else kputsn("*\t0\t0", 5, str);
   kputc('\t', str);

   // print SEQ and QUAL
   if (p.flag & 0x100) {  // for secondary alignments, don't write SEQ and QUAL
      kputsn("*\t*", 3, str);
   } else if (p.is_rev) { // the reverse strand

      // SEQ
      int i, qb = 0, qe = s->l_seq0;
      if (p.n_cigar && !is_primary && !(opt->flag&MEM_F_SOFTCLIP) && !p.is_alt) { // hard clip
         if ((p.cigar[0]&0xf) == 4 || (p.cigar[0]&0xf) == 3) qe -= p.cigar[0]>>4;
         if ((p.cigar[p.n_cigar-1]&0xf) == 4 || (p.cigar[p.n_cigar-1]&0xf) == 3) qb += p.cigar[p.n_cigar-1]>>4;
      }
      ks_resize(str, str->l + (qe - qb) + 1);
      for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq0[i]];
      kputc('\t', str);

      // QUAL
      if (s->qual) {
         ks_resize(str, str->l + (qe - qb) + 1);
         for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
         str->s[str->l] = 0;
      } else kputc('*', str);
    
   } else {                // the forward strand

      // SEQ
      int i, qb = 0, qe = s->l_seq0;
      if (p.n_cigar && !is_primary && !(opt->flag&MEM_F_SOFTCLIP) && !p.is_alt) { // hard clip
         if ((p.cigar[0]&0xf) == 4 || (p.cigar[0]&0xf) == 3) qb += p.cigar[0]>>4;
         if ((p.cigar[p.n_cigar-1]&0xf) == 4 || (p.cigar[p.n_cigar-1]&0xf) == 3) qe -= p.cigar[p.n_cigar-1]>>4;
      }
      ks_resize(str, str->l + (qe - qb) + 1);
      for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq0[i]];
      kputc('\t', str);

      // QUAL
      if (s->qual) {
         ks_resize(str, str->l + (qe - qb) + 1);
         for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
         str->s[str->l] = 0;
      } else kputc('*', str);
   }

   // TAGS
   if (p.n_cigar) {
      kputsn("\tNM:i:", 6, str); kputw(p.NM, str); // true mismatches
      // position of actual mismatches
      kputsn("\tMD:Z:", 6, str); kputs((char*)(p.cigar + p.n_cigar), str);
      kputsn("\tZC:i:", 6, str); kputw(p.ZC, str); // count of conversion
      kputsn("\tZR:i:", 6, str); kputw(p.ZR, str); // count of retention
   }
   // AS: best local SW score
   if (p.score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p.score, str); }
   // XS: 2nd best SW score or SW score of tandem hit whichever is higher
   if (p.sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(max(p.sub, p.csub), str); }
   // RG: read group
   if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
   // SA: other parts of a chimeric primary mapping
   if (regs0) mem_alnreg_tagSA(opt, bns, pac, s, p0, regs0, str);
   // PA: ratio of score / alt_score, higher the ratio, the more accurate the position
   if (is_primary && p.alt_sc > 0) ksprintf(str, "\tPA:f:%.3f", (double) p.score / p.alt_sc); // used to be lowercase pa, just to be consistent
   // XL: read length excluding adaptor
   kputsn("\tXL:i:", 6, str); kputw(s->l_seq, str);
   // XA and XB: alternative (secondary) alignment
   if (regs0) mem_alnreg_tagXAXB(opt, bns, pac, s, p0, regs0, str);
   if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
   // XR: reference/chromosome annotation
   if ((opt->flag&MEM_F_REF_HDR) && p.rid >= 0 && bns->anns[p.rid].anno != 0 && bns->anns[p.rid].anno[0] != 0) {
      int tmp;
      kputsn("\tXR:Z:", 6, str);
      tmp = str->l;
      kputs(bns->anns[p.rid].anno, str);
      unsigned i;
      for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
         if (str->s[i] == '\t') str->s[i] = ' ';
   }
   // YD: Bisulfite conversion strand label, f for forward and r for reverse, a la BWA-meth
   kputsn("\tYD:A:", 6, str);
   if (p.bss_u) kputc('u', str);
   else kputc("fr"[p.bss], str);

   kputc('\n', str);
}

/****************************************
 * output SAM format in bseq1_t *s->sam *
 ****************************************/

typedef kvec_t(int) int_v;

// select which mem_alnreg_v to output
static int_v mem_alnreg_select_format(
   const mem_opt_t *opt, const bntseq_t *bns,
   const uint8_t *pac, bseq1_t *s, mem_alnreg_v *regs) {

   // set cigar, mapq etc.
   // only the first mapping is the primary mapping
   int l; unsigned k;
   int_v to_output;  kv_init(to_output);
   for (k = l = 0; k < regs->n; ++k) {
      mem_alnreg_t *p = regs->a + k;

      // skip unmapped
      if (p->rb < 0 || p->re < 0) continue;

      // skip region lower in score than opt->T
      if (p->score < opt->T) continue;

      // skip ALT secondary mapping
      if (p->secondary >= 0 && (p->is_alt || !(opt->flag & MEM_F_ALL))) continue;

      // skip secondary mapping with score much lower than the primary mapping 
      if (p->secondary >= 0 && p->secondary < INT_MAX 
          && p->score < regs->a[p->secondary].score * opt->drop_ratio) continue;

      // flag
      // The suboptimal primary alignments are multiple (0x10000) or supplementary (0x800)
      if (l && p->secondary < 0) p->flag |= (opt->flag&MEM_F_NO_MULTI) ? 0x10000 : 0x800;
      // secondary mapping
      if (p->secondary >= 0) p->flag |= 0x100;

      // mapQ
      // set mapQ for primary mapping
      p->mapq = p->secondary < 0 ? mem_approx_mapq_se(opt, p) : 0;
      // mapq of secondary/supplementary alignment is capped by the primary mapping
      if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && l && !p->is_alt)
         p->mapq = min(p->mapq, regs->a[0].mapq);

      mem_alnreg_setSAM(opt, bns, pac, s, p);
      kv_push(int, to_output, k);
      ++l;
   }

   return to_output;
}

// Single-End
// universal_mreg is the arbitrarily selected mate reg to pair with every reg in regs
void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *regs) {

  /* if (!(opt->flag & MEM_F_ALL)) // output all alignments, hence no need to output alternatives */
  /*   mem_gen_alt(opt, bns, pac, s, regs); */
  // mem_gen_sa

  kstring_t str = {0,0,0};
  int_v to_output = mem_alnreg_select_format(opt, bns, pac, s, regs);

  // string output to s->sam
  if (to_output.n > 0) {
    // output, note each read's output depends on the cigar of other reads
    unsigned i;
    for (i = 0; i < to_output.n; ++i)
      mem_alnreg_formatSAM(opt, bns, pac, &str, s, &regs->a[to_output.a[i]], NULL, regs, !i, NULL);
  } else { // unmapped read
    mem_alnreg_t reg; memset(&reg, 0, sizeof(mem_alnreg_t));
    reg.rid = -1;
    reg.flag = 0x4;
    mem_alnreg_formatSAM(opt, bns, pac, &str, s, &reg, NULL, regs, 1, NULL);
  }
  s->sam = str.s;
  kv_destroy(to_output);
}

// Paired-End
// pairing with the best in mate alignment
void mem_reg2sam_pe_nopairing(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t s[2], mem_alnreg_v regs_pair[2], mem_pestat_t pes) {

  mem_alnreg_t *best[2] = {0,0};
  int_v to_outputs[2];
  mem_alnreg_t unmapped_reg[2]; // the "unmapped" alignment

  if (bwa_verbose >= 4) printf("PE no pairing.\n");

  // looking for the best alnreg to pair
  int i;
  for (i = 0; i < 2; ++i) {
    mem_alnreg_v *regs = regs_pair + i;
    to_outputs[i] = mem_alnreg_select_format(opt, bns, pac, &s[i], regs);
    if (to_outputs[i].n > 0) {
      best[i] = &regs->a[to_outputs[i].a[0]];
    } else {
      memset(&unmapped_reg[i], 0, sizeof(mem_alnreg_t));
      unmapped_reg[i].rid = -1; 
      unmapped_reg[i].flag = 0x40<<i | 0x1 | 0x4;
      best[i] = &unmapped_reg[i];
    }
  }

  // output
  for (i = 0; i < 2; ++i) {
    kstring_t str = {0,0,0};
    mem_alnreg_v *regs = regs_pair + i;
    if (to_outputs[i].n) {
      unsigned j;
      for (j = 0; j < to_outputs[i].n; ++j) {
        mem_alnreg_t *p = &regs->a[to_outputs[i].a[j]];
        if (!best[!i]) p->flag |= 0x8; // distinguish unmapped mate from unpaired read
        mem_alnreg_formatSAM(opt, bns, pac, &str, &s[i], &regs->a[to_outputs[i].a[j]], best[!i], regs, !j, &pes);
      }
    } else mem_alnreg_formatSAM(opt, bns, pac, &str, &s[i], best[i], best[!i], NULL, 1, &pes);

    s[i].sam = str.s;
  }

  for (i = 0; i < 2; ++i) kv_destroy(to_outputs[i]);
}

#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))
void mem_reg2sam_pe(
   const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
   uint64_t id, bseq1_t s[2], mem_alnreg_v regs_pair[2], mem_pestat_t pes) {

   if (bwa_verbose >= 4) {
      printf("[%s] Read 1 in pairing:\n", __func__);
      mem_print_regions(bns, &regs_pair[0]);
      printf("[%s] Read 2 in pairing:\n", __func__);
      mem_print_regions(bns, &regs_pair[1]);
      printf("\n");
   }

   // flags for paired reads
   int i; unsigned k;
   for (i=0; i<2; ++i)
      for (k = 0; k<regs_pair[i].n; ++k)
         regs_pair[i].a[k].flag |= (0x40 << i) | 1; // set which is read1/2
  
   if (opt->flag & MEM_F_NOPAIRING)
      return mem_reg2sam_pe_nopairing(opt, bns, pac, s, regs_pair, pes);
   if (regs_pair[0].n_pri == 0 || regs_pair[1].n_pri == 0)
      return mem_reg2sam_pe_nopairing(opt, bns, pac, s, regs_pair, pes);

   // check if an end has multiple hits even after mate-SW, skip pairing if the case
   int is_multi[2]; unsigned j;
   for (i = 0; i < 2; ++i) {
      for (j = 1; j < regs_pair[i].n_pri; ++j) {// start from the second
         //printf("final %d:%u - %d - %d, %d\n",i,j,regs_pair[i].a[j].secondary,regs_pair[i].a[j].score, regs_pair[i].n_pri); 
         //mem_print_region1(bns, &regs_pair[i].a[j]);
         if (regs_pair[i].a[j].secondary < 0 && regs_pair[i].a[j].score >= opt->T) break;
      }
      // if there is a primary chromosome, primary, good alignment
      is_multi[i] = j < regs_pair[i].n_pri ? 1 : 0;
   }
   // TODO: in rare cases, the true hit may be long but with low score
   if (is_multi[0] || is_multi[1]) return mem_reg2sam_pe_nopairing(opt, bns, pac, s, regs_pair, pes);

   // Actual pairing and set mapQ
   int pscore, sub_pscore; // best and 2nd best pairing score
   int n_subpairings; int z[2] = {0};
   mem_pair(opt, bns, pes, regs_pair, id, &pscore, &sub_pscore, &n_subpairings, z);
   if (pscore <= 0) return mem_reg2sam_pe_nopairing(opt, bns, pac, s, regs_pair, pes);

   if (bwa_verbose >= 4) {
      mem_alnreg_t *p1 = &regs_pair[0].a[z[0]];
      mem_alnreg_t *p2 = &regs_pair[1].a[z[1]];
      mem_alnreg_setSAM(opt, bns, pac, &s[0], p1);
      mem_alnreg_setSAM(opt, bns, pac, &s[1], p2);
      printf("** pairing read 1: %d, [%d,%d) <=> [%ld,%ld,%s,%d) <> read 2: %d, [%d,%d) <=> [%ld,%ld,%s,%d)\n", 
             p1->score, p1->qb, p1->qe, (long)p1->rb, (long)p1->re, bns->anns[p1->rid].name, p1->pos,
             p2->score, p2->qb, p2->qe, (long)p2->rb, (long)p2->re, bns->anns[p2->rid].name, p2->pos);
   }
  
   // opt->pen_unpaired - penalty for not pairing
   int score_unpaired = regs_pair[0].a[0].score + regs_pair[1].a[0].score - opt->pen_unpaired;
   if (pscore > score_unpaired) { // use z for pairing
      // mapQ of pairing (q_pe)
      if (bwa_verbose >= 4) printf("Favor pairing\n");
      sub_pscore = max(sub_pscore, score_unpaired);
      int q_pe = raw_mapq(pscore - sub_pscore, opt->a);
      if (n_subpairings > 0) q_pe -= (int)(4.343 * log(n_subpairings+1) + .499);
      q_pe = max(0, min(60, q_pe));
      q_pe = (int)(q_pe * (1. - .5 * (regs_pair[0].a[0].frac_rep + regs_pair[1].a[0].frac_rep)) + .499);

      // mapQ of each read when paired
      int q_se[2]; // mapQ of each read in the pair
      mem_alnreg_t *c[2];
      c[0] = &regs_pair[0].a[z[0]];
      c[1] = &regs_pair[1].a[z[1]];
      for (i = 0; i < 2; ++i) {
         if (c[i]->secondary >= 0) {
            c[i]->sub = regs_pair[i].a[c[i]->secondary].score;
            c[i]->secondary = -2;
         }
         q_se[i] = mem_approx_mapq_se(opt, c[i]);
      }

      // add 40 or reach q_pe, whichever is smaller
      q_se[0] = max(q_se[0], min(q_pe, q_se[0] + 40));
      q_se[1] = max(q_se[1], min(q_pe, q_se[1] + 40));

      // cap at the tandem repeat score
      c[0]->mapq = min(q_se[0], raw_mapq(c[0]->score - c[0]->csub, opt->a));
      c[1]->mapq = min(q_se[1], raw_mapq(c[1]->score - c[1]->csub, opt->a));
   } else {                         // use best hits for pairing
      if (bwa_verbose >= 4) printf("Favor best hits in pairing\n");
      z[0] = z[1] = 0;
      regs_pair[0].a[0].mapq = mem_approx_mapq_se(opt, &regs_pair[0].a[0]);
      regs_pair[1].a[0].mapq = mem_approx_mapq_se(opt, &regs_pair[1].a[0]);
   }

   // if the chosen read is a secondary, switch it with its designated primary
   for (i = 0; i < 2; ++i) {
      mem_alnreg_v *regs = &regs_pair[i];
      int k = regs->a[z[i]].secondary_all;
      // switch secondary and primary if both of them are non-ALT
      if (k >= 0 && (unsigned) k < regs->n_pri) {
         // the old primary is not a secondnary by itself
         assert(regs->a[k].secondary_all < 0);
         for (j = 0; j < regs->n; ++j) {
            // make the old primary and every secondary of
            // the old primary the secondnary of the chosen
            if (regs->a[j].secondary_all == k || j == (unsigned) k)
               regs->a[j].secondary_all = z[i];
         }
         regs->a[z[i]].secondary_all = -1; // the chosen is now a primary
      }
   }

   // set SAM
   for (i = 0; i < 2; ++i) {
      mem_alnreg_setSAM(opt, bns, pac, &s[i], &regs_pair[i].a[z[i]]);
   }

   // write SAM
   for (i = 0; i < 2; ++i) {
      kstring_t str = {0,0,0};
      mem_alnreg_v *regs = &regs_pair[i];
      mem_alnreg_t *reg = regs_pair[i].a + z[i];
      mem_alnreg_t *mreg = regs_pair[!i].a + z[!i];

      mem_alnreg_formatSAM(opt, bns, pac, &str, &s[i], reg, mreg, regs, 1, &pes);

      // output one ALT hit as unpaired mapping?
      if (regs->n_pri < regs->n) {
         mem_alnreg_t *p = &regs->a[regs->n_pri]; // output the best ALT hit
         if (p->score >= opt->T && p->secondary < 0) {
            p->flag |= 0x800; // supplementary alignment
            mem_alnreg_setSAM(opt, bns, pac, &s[i], p);
            mem_alnreg_formatSAM(opt, bns, pac, &str, &s[i], p, NULL, regs, 0, &pes); // is mate none?
         }
      }
      s[i].sam = str.s;
   }
}
