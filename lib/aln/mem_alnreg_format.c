/* SAM formating for mem_regaln_t
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
#include <stdlib.h>
#include "bwamem.h"
#include "kvec.h"
#include "kstring.h"
#include "wzmisc.h"
#include "sam.h"


// set CIGAR
// mapQ and flag are set outside this function
static void mem_alnreg_setSAM(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_t *reg) {

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

  // incrementally double bandwidth
  uint32_t *cigar = 0; int n_cigar;
  int score; int last_sc = -(1<<30);
  for (i=0; i<3; ++i, w<<=1, last_sc=score) {
    free(cigar);
    w = min(w, opt->w<<2);

    // regenerate cigar related info with new bandwidth
    cigar = bis_bwa_gen_cigar2(reg->parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, pac, reg->qe - reg->qb, (uint8_t*) &query[reg->qb], reg->rb, reg->re, &score, &n_cigar, &reg->NM, &reg->ZC, &reg->ZR, reg->parent);

    if (bwa_verbose >= 4) printf("* Final alignment: w=%d, global_sc=%d, local_sc=%d\n", w, score, reg->truesc);

    if (score == last_sc) break;
    if (w == opt->w << 2) break;
    if (score >= reg->truesc - opt->a) break;
  }
  int l_MD = strlen((char*) (cigar+n_cigar))+1;

  // pos
  int is_rev;
  int64_t rpos = bns_depos(bns, reg->rb < bns->l_pac ? reg->rb : reg->re-1, &is_rev);

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
  if (reg->qb != 0 || reg->qe != s->l_seq) {
    int clip5, clip3;
    clip5 = is_rev ? s->l_seq - reg->qe : reg->qb;
    clip3 = is_rev ? reg->qb : s->l_seq - reg->qe;
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
  reg->cigar = cigar;
  reg->n_cigar = n_cigar;

  assert(bns_pos2rid(bns, rpos) == reg->rid);
  reg->pos = rpos - bns->anns[reg->rid].offset;
  reg->sam_set = 1;

  return;
}

static void mem_alnreg_freeSAM(mem_alnreg_v *regs) {
  int i; unsigned j;
  for (i = 0; i < 2; ++i)
    for (j = 0; j < regs->n; ++j)
      free(regs->a[j].cigar);
}

/* Generate XA, put to str */
static void mem_alnreg_setXA(const mem_opt_t *opt, const bntseq_t *bns, const mem_alnreg_t *p0, const mem_alnreg_v *regs0, kstring_t *sam_str) {

  if (!regs0 || !(opt->flag & MEM_F_ALL)) return;

  int cnt = 0, has_alt = 0; unsigned i;
  for (i=0; i<regs0->n; ++i) {
    int r = get_pri_idx(opt->XA_drop_ratio, regs0->a, i);
    if (r >= 0 && regs0->a + r == p0) {
      ++cnt;
      if (regs0->a[i].is_alt) has_alt = 1;
    }
  }
  
  // skip if the total number of alts is higher than opt->max_XA_hits_alt
  // or if all primary number of alts is higher than opt->max_XA_hits
  if (cnt > opt->max_XA_hits_alt || (!has_alt && cnt > opt->max_XA_hits)) return;

  kstring_t str = {0};
  for (i=0; i<regs0->n; ++i) {

    mem_alnreg_t *q = regs0->a + i;
    int r = get_pri_idx(opt->XA_drop_ratio, regs0->a, i);
    if (r < 0 || regs0->a + r != p0 || q->n_cigar == 0) continue;

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
    kputc(';', &str);
  }

  if (str.l) {
    kputsn("\tXA:Z:", 6, sam_str); 
    kputs(str.s, sam_str);
  }
  free(str.s);
}

/* Generate SA-tag, put to str */
static void mem_alnreg_setSA(const bntseq_t *bns, const mem_alnreg_t *p0, const mem_alnreg_v *regs0, kstring_t *sam_str) {

  if (!regs0 || p0->flag & 0x100) return;

  kstring_t str = {0};
  unsigned i;
  for (i=0; i < regs0->n; ++i) {
    mem_alnreg_t *q = regs0->a + i;
    if (q == p0 || q->n_cigar == 0 || q->flag & 0x100) continue;

    kputs(bns->anns[q->rid].name, &str);
    kputc(',', &str);
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

int is_proper_pair(mem_alnreg_t *r1, mem_alnreg_t *r2, mem_pestat_t *pes) {

  // switch 1 and 2 if flag indicates otherwise
  if (r1->flag & BAM_FREAD2 && r2->flag & BAM_FREAD1) {
    mem_alnreg_t *tmp = r2;
    r2 = r1; r1 = tmp;
  }
    
  if (r1->rid != r2->rid) return 0;
  if (r1->pos - r2->pos < pes->high && r1->pos - r2->pos > pes->low) return 1;
  else return 0;
}

/*************************
 * format SAM
 *************************
 * mate is set at final stage because the mate might be asymmetric with alternative mappings
 * It doesn't not change the mem_alnreg_t inputs
 *************************/
void mem_alnreg_formatSAM(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, const mem_alnreg_t *p0, const mem_alnreg_t *m0, const mem_alnreg_v *regs0, int is_primary, mem_pestat_t *pes) {

  // make copies
  mem_alnreg_t p = *p0;
  mem_alnreg_t m = {0};
  if (m0) m = *m0;

  // set mate-related flags
  p.flag |= m0 ? 0x1 : 0; // is paired in sequencing
  p.flag |= m0 && m.rid < 0 ? 0x8 : 0; // is mate mapped

  if (p.rid >= 0 && m0 && m.rid >= 0 && pes && is_proper_pair(&p, &m, pes)) {
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
  ks_resize(str, str->l + s->l_seq + l_name + (s->qual ? s->l_seq : 0) + 20);
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
      int64_t p0 = p.pos + (p.is_rev? get_rlen(p.n_cigar, p.cigar) - 1 : 0);
      int64_t p1 = m.pos + (m.is_rev? get_rlen(m.n_cigar, m.cigar) - 1 : 0);
      if (m.n_cigar == 0 || p.n_cigar == 0) kputc('0', str);
      else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
    } else kputc('0', str);
  } else kputsn("*\t0\t0", 5, str);
  kputc('\t', str);

  // print SEQ and QUAL
  if (p.flag & 0x100) {  // for secondary alignments, don't write SEQ and QUAL
    kputsn("*\t*", 3, str);
  } else if (p.is_rev) { // the reverse strand

    // SEQ
    int i, qb = 0, qe = s->l_seq;
    if (p.n_cigar && !is_primary && !(opt->flag&MEM_F_SOFTCLIP) && !p.is_alt) { // hard clip
      if ((p.cigar[0]&0xf) == 4 || (p.cigar[0]&0xf) == 3) qe -= p.cigar[0]>>4;
      if ((p.cigar[p.n_cigar-1]&0xf) == 4 || (p.cigar[p.n_cigar-1]&0xf) == 3) qb += p.cigar[p.n_cigar-1]>>4;
    }
    ks_resize(str, str->l + (qe - qb) + 1);
    for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
    kputc('\t', str);

    // QUAL
    if (s->qual) {
      ks_resize(str, str->l + (qe - qb) + 1);
      for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
      str->s[str->l] = 0;
    } else kputc('*', str);
    
  } else {                // the forward strand

    // SEQ
    int i, qb = 0, qe = s->l_seq;
    if (p.n_cigar && !is_primary && !(opt->flag&MEM_F_SOFTCLIP) && !p.is_alt) { // hard clip
      if ((p.cigar[0]&0xf) == 4 || (p.cigar[0]&0xf) == 3) qb += p.cigar[0]>>4;
      if ((p.cigar[p.n_cigar-1]&0xf) == 4 || (p.cigar[p.n_cigar-1]&0xf) == 3) qe -= p.cigar[p.n_cigar-1]>>4;
    }
    ks_resize(str, str->l + (qe - qb) + 1);
    for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
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
    kputsn("\tNM:i:", 6, str); kputw(p.NM, str);
    kputsn("\tMD:Z:", 6, str); kputs((char*)(p.cigar + p.n_cigar), str);
    kputsn("\tZC:i:", 6, str); kputw(p.ZC, str);
    kputsn("\tZR:i:", 6, str); kputw(p.ZR, str);
  }
  // AS: best local SW score
  if (p.score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p.score, str); }
  // XS: 2nd best SW score or SW score of tandem hit whichever is higher
  if (p.sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(max(p.sub, p.csub), str); }
  // RG: read group
  if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
  // SA: other parts of a chimeric primary mapping
  if (regs0) mem_alnreg_setSA(bns, p0, regs0, str);
  // PA: ratio of score / alt_score, higher the ratio, the more accurate the position
  if (is_primary) ksprintf(str, "\tPA:f:%.3f", (double) p.score / p.alt_sc); // used to be lowercase pa, just to be consistent
  // XA: alternative alignment
  if (regs0) mem_alnreg_setXA(opt, bns, p0, regs0, str);
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
  // YD: Bisulfite conversion strand label, per BWA-meth
  kputsn("\tYD:A:", 6, str);
  if (p.bss < 0) kputc('u', str);
  else kputc("fr"[p.bss], str);

  kputc('\n', str);
}

/****************************************
 * output SAM format in bseq1_t *s->sam *
 ****************************************/

// Single-End
// universal_mreg is the arbitrarily selected mate reg to pair with every reg in regs
void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *regs, const mem_alnreg_t *universal_mreg) {

  if (!(opt->flag & MEM_F_ALL)) // output all alignments, hence no need to output alternatives
    mem_gen_alt(opt, bns, pac, s, regs);
  // mem_gen_sa

  kstring_t str = {0};

  // set cigar, mapq etc.
  // only the first mapping is the primary mapping
  int l; unsigned k;
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
    // is on the reverse strand
    p->flag |= p->is_rev ? 0x10 : 0;

    // mapQ
    // set mapQ for primary mapping
    p->mapq = p->secondary < 0 ? mem_approx_mapq_se(opt, p) : 0;
    // mapq of secondary/supplementary alignment is capped by the primary mapping
    if (l && !p->is_alt) p->mapq = min(p->mapq, regs->a[0].mapq);

    mem_alnreg_setSAM(opt, bns, pac, s, p);
    mem_alnreg_formatSAM(opt, bns, &str, s, p, universal_mreg, regs, !k, NULL);

    ++l;
  }
  mem_alnreg_freeSAM(regs);

  // string output to s->sam
  if (l == 0) { // unmapped read
    mem_alnreg_t reg = {0};
    reg.rid = -1;
    reg.flag |= 0x4;
    mem_alnreg_formatSAM(opt, bns, &str, s, &reg, universal_mreg, regs, 1, NULL);
  }
  s->sam = str.s;
}

// Paired-End
#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))
void mem_reg2sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint64_t id, bseq1_t s[2], mem_alnreg_v regs_pair[2], mem_pestat_t *pes) {

  // flags for paired reads
  int i; unsigned k;
  for (i=0; i<2; ++i)
    for (k = 0; k<regs_pair[i].n; ++k)
      regs_pair[i].a[k].flag |= (0x40 << i) | 1; // set which is read1/2
  
  if (opt->flag & MEM_F_NOPAIRING) goto NO_PAIRING;
  if (regs_pair[0].n_pri == 0 || regs_pair[1].n_pri == 0) goto NO_PAIRING;

  // check if an end has multiple hits even after mate-SW, skip pairing if the case
  int is_multi[2]; unsigned j;
  for (i = 0; i < 2; ++i) {
    for (j = 1; j < regs_pair[i].n_pri; ++j) // start from the second
      if (regs_pair[i].a[j].secondary < 0 && regs_pair[i].a[j].score >= opt->T) break;
    // if there is a primary chromosome, primary, good alignment
    is_multi[i] = j < regs_pair[i].n_pri ? 1 : 0;
  }
  // TODO: in rare cases, the true hit may be long but with low score
  if (is_multi[0] || is_multi[1]) goto NO_PAIRING;

  // Actual pairing and set mapQ
  int pscore, sub_pscore; // best and 2nd best pairing score
  int n_subpairings; int z[2];
  mem_pair(opt, bns, pac, pes, s, regs_pair, id, &pscore, &sub_pscore, &n_subpairings, z);
  if (pscore <= 0) goto NO_PAIRING;
  // opt->pen_unpaired - penalty for not pairing
  int score_unpaired = regs_pair[0].a[0].score + regs_pair[1].a[0].score - opt->pen_unpaired;
  kstring_t str; str.l = str.m = 0; str.s = 0;
  if (pscore > score_unpaired) { // use z for pairing
    // mapQ of pairing (q_pe)
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
    z[0] = z[1] = 0;
    regs_pair[0].a[0].mapq = mem_approx_mapq_se(opt, &regs_pair[0].a[0]);
    regs_pair[1].a[0].mapq = mem_approx_mapq_se(opt, &regs_pair[1].a[0]);
  }

  // if the chosen read is a secondary, switch it with its designated primary
  for (i = 0; i < 2; ++i) {
    mem_alnreg_v *regs = &regs_pair[i];
    int k = regs->a[z[i]].secondary_all;
    if (k >= 0 && (unsigned) k < regs->n_pri) { // switch secondary and primary if both of them are non-ALT
      assert(regs->a[k].secondary_all < 0); // the old primary is not a secondnary by itself
      for (j = 0; j < regs->n; ++j)
        if (regs->a[j].secondary_all == k || j == (unsigned) k) // make the old primary and every secondary of the old primary the secondnary of the chosen
          regs->a[j].secondary_all = z[i];
      regs->a[z[i]].secondary_all = -1; // the chosen is now a primary
    }
  }

  // write SAM
  for (i = 0; i < 2; ++i) {

    if (!(opt->flag & MEM_F_ALL)) mem_gen_alt(opt, bns, pac, &s[i], &regs_pair[i]);

    mem_alnreg_t *reg = regs_pair[i].a + z[i];
    mem_alnreg_t *mreg = regs_pair[!i].a + z[!i];
    mem_alnreg_v *regs = &regs_pair[i];

    mem_alnreg_setSAM(opt, bns, pac, &s[i], reg);
    mem_alnreg_formatSAM(opt, bns, &str, &s[i], reg, mreg, regs, 1, pes);

    // h[i].XA = XA[i]? XA[i][z[i]] : 0;
    // aa[i][n_aa[i]++] = h[i];

    // output one ALT hit as unpaired mapping?
    if (regs->n_pri < regs->n) {
      mem_alnreg_t *p = &regs->a[regs->n_pri]; // output the best ALT hit
      if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
      p->flag |= 0x800;
      mem_alnreg_setSAM(opt, bns, pac, &s[i], p);
      mem_alnreg_formatSAM(opt, bns, &str, &s[i], p, NULL, regs, 0, pes); // is mate none?
      // g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p);
      // g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
      // aa[i][n_aa[i]++] = g[i];
    }
    s[i].sam = strdup(str.s); str.l = 0;

    mem_alnreg_freeSAM(&regs_pair[i]);
  }

  mem_alnreg_t *best[2] = {0};
 NO_PAIRING:                // pairing with the best in mate alignment

  // looking for the best alnreg to pair in the following order
  // 1) best primary chromosome alnrneg if the score > opt->T; 
  // 2) best non-primary chromosome alnreg if the score > opt->T
  // 3) best alnreg primary or non-primary
  for (i = 0; i < 2; ++i) {
    mem_alnreg_v *regs = regs_pair + i;
    best[i] = &regs->a[0];
    if (best[i]->score < opt->T && 
        regs_pair[i].n_pri < regs->n && 
        regs->a[regs_pair[i].n_pri].score >= opt->T)
      best[i] = regs->a + regs_pair[i].n_pri;
  }

  // output
  for (i = 0; i < 2; ++i)
    mem_reg2sam_se(opt, bns, pac, &s[i], &regs_pair[i], best[!i]);

  /* if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name); */
  return;
}

