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

// note: mapQ is set outside this function
void mem_alnreg_setSAM(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_t *reg, int extra_flag) {

  if (reg == 0 || reg->rb < 0 || reg->re < 0) { // unmapped record
    unmapped_read_format(sam_str);
    return;
  }

  // nt4 encoding
  uint8_t *query;
  int i;
  for (i = 0; i < s->l_seq; ++i) {
    query[i] = query_[i] < 5 ? query_[i] : nst_nt4_table[(int)query_[i]];
  }

  /** flag **/
  reg->flag = extra_flag;
  if (reg->secondary >= 0) flag |= 0x100; // secondary mapping

  /** cigar and pos **/
  // initial bandwidth is the max of insertion and deletion calculation
  w = max(
      infer_bw(reg->qe-reg->qb, reg->re-reg->rb, opt->a, opt->o_del, opt->e_del),
      infer_bw(reg->qe-reg->qb, reg->re-reg->rb, opt->a, opt->o_ins, opt->e_ins));
  if (w > opt->w) w = min(w, reg->w);

  // incrementally double bandwidth
  uint32_t *cigar = 0; int n_cigar;
  int score; int last_sc = -(1<<30);
  for (i=0; i<3; ++i, w<<=1, last_sc=score) {
    free(cigar);
    w = min(w, opt->w<<2);

    // regenerate cigar related info with new bandwidth
    cigar = bis_bwa_gen_cigar2(reg->parent?opt->ctmat:opt->gamat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, reg->qe - reg->qb, (uint8_t*) &query[qb], reg->rb, reg->re, &score, &n_cigar, &reg->NM, &reg->ZC, &reg->ZR, reg->parent);

    if (bwa_verbose >= 4) printf("* Final alignment: w=%d, global_sc=%d, local_sc=%d\n", w, score, reg->truesc);

    if (score = last_sc) break;
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
  if (qb != 0 || qe != l_query) {
    int clip5, clip3;
    clip5 = is_rev ? l_query - qe : qb;
    clip3 = is_rev ? qb : l_query - qe;
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
  reg->pos = rpos - bns->anns[rid].offset;
  reg->sam_set = 1;
  return;
}

void mem_alnreg_formatSAM(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, const mem_alnreg_t *p, const mem_alnreg_t *m, int is_primary) {

  if (!reg->sam_set) return;

  // set flag
  p->flag |= m ? 0x1 : 0; // is paired in sequencing
  p->flag |= p->rid < 0 ? 0x4 : 0; // is mapped, should always be true here?
  p->flag |= m && m->rid < 0 ? 0x8 : 0; // is mate mapped

  // copy mate coordinate to alignment
  if (p->rid < 0 && m && m->rid >= 0) {
    p->rid = m->rid;
    p->pos = m->pos;
    p->is_rev = m->is_rev; // not 100% sure if we should copy is_rev
    p->n_cigar = 0;
  }

  // copy alignment coordinates to mate
  int m_rid = m->rid;
  int m_pos = m->pos;
  int m_n_cigar = m->n_cigar;
  unsigned m_is_rev = m->is_rev;
  if (m && (!m->sam_set && m->rid < 0) && p->rid >= 0) {  // copy alignment to mate
    m_rid = p->rid;
    m_pos = p->pos;
    m_is_rev = p->is_rev;
    m_n_cigar = 0;
  }
  p->flag |= p->is_rev ? 0x10 : 0; // is on the reverse strand
  p->flag |= m && m_is_rev ? 0x20 : 0; // is mate on the reverse strand

  // print up to CIGAR
  int l_name = strlen(s->name);
  ks_resize(str, str->l + s->l_seq + l_name + (s->qual ? s->l_seq : 0) + 20);
  kputsn(s->name, l_name, str); kputc('\t', str); // read name, QNAME
  kputw((p->flag & 0xffff) | (p->flag & 0x10000 ? 0x100 : 0), str); kputc('\t', str); // FLAG
  if (p->rid >= 0) { // with coordinate
    kputs(bns->anns[p->rid].name, str); kputc('\t', str); // reference/chromosome name, RNAME
    kputl(p->pos + 1, str); kputc('\t', str); // POS
    kputw(p->mapq, str); kputc('\t', str); // MAPQ
    if (p->n_cigar) { // CIGAR
      for (i = 0; i < p->n_cigar; ++i) {
        int c = p->cigar[i] & 0xf;
        if (!(opt->flag & MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
          c = is_primary ? 3 : 4; // use hard clipping for supplementary alignments
        kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
      }
    } else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
  } else kputsn("*\t0\t0\t*", 7, str); // without coordinte
  kputc('\t', str);

  // print the mate position if applicable
  if (m && m_rid >= 0) {
    if (p->rid == m_rid) kputc('=', str);
    else kputs(bns->anns[m_rid].name, str);
    kputc('\t', str);
    kputl(m_pos + 1, str); kputc('\t', str);
    if (p->rid == m_rid) {
      int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
      int64_t p1 = m_pos + (m_is_rev? get_rlen(m_n_cigar, m_cigar) - 1 : 0);
      if (m_n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
      else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
    } else kputc('0', str);
  } else kputsn("*\t0\t0", 5, str);
  kputc('\t', str);

  // print SEQ and QUAL
  if (p->flag & 0x100) {  // for secondary alignments, don't write SEQ and QUAL
    kputsn("*\t*", 3, str);
  } else if (p->is_rev) { // the reverse strand

    // SEQ
    int i, qb = 0, qe = s->l_seq;
    if (p->n_cigar && !is_primary && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // hard clip
      if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
      if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
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
    if (p->n_cigar && !is_primary && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // hard clip
      if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
      if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
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
  if (p->n_cigar) {
    kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
    kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
    kputsn("\tZC:i:", 6, str); kputw(p->ZC, str);
    kputsn("\tZR:i:", 6, str); kputw(p->ZR, str);
  }
  // AS: best local SW score
  if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
  // XS: 2nd best SW score or SW score of tandem hit whichever is higher
  if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(max(p->sub, p->csub), str); }
  // RG: read group
  if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
  // SA: other primary hits
  if (is_primary && p->SA) { kputsn("\tSA:Z:", 6 , str); kputs(p->SA, str);}
  // PA: ratio of score / alt_score, higher the ratio, the more accurate the position
  if (is_primary) ksprintf(str, "\tPA:f:%.3f", (double) p->score / p->alt_sc); // used to be lowercase pa, just to be consistent
  // XA: alternative alignment
  if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
  if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
  // XR: reference/chromosome annotation
  if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
    int tmp;
    kputsn("\tXR:Z:", 6, str);
    tmp = str->l;
    kputs(bns->anns[p->rid].anno, str);
    for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
      if (str->s[i] == '\t') str->s[i] = ' ';
  }
  // YD: Bisulfite conversion strand label, per BWA-meth
  kputsn("\tYD:A:", 6, str);
  if (p->bss < 0) kputc('u', str);
  else kputc("fr"[p->bss], str);

  kputc('\n', str);
}

/****************************************
 * output SAM format in bseq1_t *s->sam *
 ****************************************/

// Single-End
void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *regs, int extra_flag) {
  
  if (!(opt->flag & MEM_F_ALL)) // output all alignments, hence no need to output alternatives
    XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
  // mem_gen_sa

  // set cigar, mapq etc.
  // only the first mapping is the primary mapping
  unsigned k; int l;
  for (k = l = 0; k < regs->n; ++k) {
    mem_alnreg_t *p = regs->a + k;

    // skip region lower in score than opt->T
    if (p->score < opt->T) continue;

    // skip non-primary chromosome secondary mapping
    if (p->secondary >= 0 && (p->is_alt || !(opt->flag & MEM_F_ALL))) continue;

    // skip secondary mapping with score much lower than the primary mapping 
    if (p->secondary >= 0 && p->secondary < INT_MAX 
        && p->score < regs->a[p->secondary].score * opt->drop_ratio) continue;

    // set mapQ
    p->mapq = p->secondary < 0 ? mem_approx_mapq_se(opt, p) : 0;

    mem_alnreg_setSAM(opt, bns, pac, s, p, extra_flag);

    // keep only 1 primary alignment, others are either secondary or supplementary
    if (l && p->secondary < 0) q->flag |= (opt->flag&MEM_F_NO_MULTI) ? 0x10000 : 0x800;
    
    // mapq of secondary/supplementary alignment is capped by the primary mapping
    if (l && !p->is_alt) p->mapq = min(p->mapq, regs->a[0].mapq);

    ++l;
  }

  // string output to s->sam
  kstring_t sam_str;
  sam_str.l = sam_str.m = 0; sam_str.s = 0;
  if (l == 0) { // unmapped read
    mem_alnreg_t reg = {0}; reg.rid = -1; reg.sam_set = 1;
    mem_alnreg_formatSAM(opt, bns, &str, s, &reg, NULL, 1);
  } else {
    for (k = 0; k < regs->n; ++k) {
      mem_alnreg_formatSAM(opt, bns, &str, s, &regs->a[k], NULL, !k);
      free(regs->a[k].cigar);
    }
  }
  s->sam = str.s;
}

// Paired-End
void mem_reg2sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint64_t id, bseq1_t s[2], mem_alnreg_v regs_pair[2]) {

  int n = 0, i, j; 

  kstring_t str;
  str.l = str.m = 0; str.s = 0;

  if (!(opt->flag & MEM_F_NO_RESCUE)) mem_alnreg_matesw(opt, bns, pac, pes, s[2], regs_pair[2]);

  int n_pri[2];
  n_pri[0] = mem_mark_primary_se(opt, regs_pair[0].n, regs_pair[0].a, id<<1|0);
  n_pri[1] = mem_mark_primary_se(opt, regs_pair[1].n, regs_pair[1].a, id<<1|1);

  if (opt->flag & MEM_F_NOPAIRING) goto NO_PAIRING;

  int extra_flag = 1;

  /* pairing mate reads */
  int pscore, sub_pscore; // best and 2nd best pairing score
  int n_subpairings; int z[2];
  mem_pair(opt, bns, pac, pes, s, regs_pair, id, &pscore, &sub_pscore, &n_subpairing, z, n_pri);
  if (n_pri[0] && n_pri[1] && pscore > 0) {

    // check if an end has multiple hits even after mate-SW
    int is_multi[2];
    for (i = 0; i < 2; ++i) {
      for (j = 1; j < n_pri[i]; ++j)
        if (regs_pair[i].a[j].secondary < 0 && regs_pair[i].a[j].score >= opt->T) break;
      is_multi[i] = j < n_pri[i]? 1 : 0;
    }
    if (is_multi[0] || is_multi[1]) goto NO_PAIRING; // TODO: in rare cases, the true hit may be long but with low score

    // The following decides whether to pair.
    // It also sets mapQ assuming no split hits.
    int score_unpaired = regs_pair[0].a[0].score + regs_pair[1].a[0].score - opt->pen_unpaired; // opt->pen_unpaired - penalty for not pairing
    if (pscore > score_unpaired) { // paired alignment is preferred

      // mapQ of pairing
      sub_pscore = max(sub_pscore, score_unpaired);
      int q_pe = raw_mapq(pscore - sub_pscore, opt->a); // mapQ for pairing
      if (n_subpairings > 0) q_pe -= (int)(4.343 * log(n_subpairings+1) + .499);
      if (q_pe < 0) q_pe = 0;
      if (q_pe > 60) q_pe = 60;
      q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499);

      // mapQ of each read when paired
      int q_se[2]; // mapQ of each read in the pair
      mem_alnreg_t *c[2];
      c[0] = &regs_pair[0].a[z[0]]; c[1] = &regs_pair[1].a[z[1]];
      for (i = 0; i < 2; ++i) {
        if (c[i]->secondary >= 0) {
          c[i]->sub = regs_pair[i].a[c[i]->secondary].score;
          c[i]->secondary = -2;
        }
        q_se[i] = mem_approx_mapq_se(opt, c[i]);
      }
      q_se[0] = q_se[0] > q_pe? q_se[0] : q_pe < q_se[0] + 40? q_pe : q_se[0] + 40;
      q_se[1] = q_se[1] > q_pe? q_se[1] : q_pe < q_se[1] + 40? q_pe : q_se[1] + 40;
      extra_flag |= 2;
      // cap at the tandem repeat score
      q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a) ? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
      q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a) ? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);
      c[0]->mapq = q_se[0];
      c[1]->mapq = q_se[1];
    } else {                         // the unpaired alignment is preferred
      z[0] = z[1] = 0;
      regs_pair[0].a[0].mapq = mem_approx_mapq_se(opt, &regs_pair[0].a[0]);
      regs_pair[1].a[0].mapq = mem_approx_mapq_se(opt, &regs_pair[1].a[0]);
    }

    // if the chosen read is a secondary, switch it with its designated primary
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

    char **XA[2];
    if (!(opt->flag & MEM_F_ALL)) {
      for (i = 0; i < 2; ++i)
        XA[i] = mem_gen_alt(opt, bns, pac, &a[i], s[i].l_seq, s[i].seq);
    } else XA[0] = XA[1] = 0;

    /* write SAM */
    for (i = 0; i < 2; ++i) {
      mem_alnreg_setSAM(opt, bns, pac, &s[i], &regs_pair[i].a[z[i]], 0x40<<i | extra_flag);
      mem_alnreg_formatSAM(opt, bns, &str, &s[i], &regs_pair[i].a[z[i]], &regs_pair[i].a[z[!i]], 1);
      // h[i].XA = XA[i]? XA[i][z[i]] : 0;
      // aa[i][n_aa[i]++] = h[i];
      if (n_pri[i] < regs_pair[i].n) { // the read has ALT hits
        mem_alnreg_t *p = &regs_pair[i].a[n_pri[i]]; // output the best ALT hit
        if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
        mem_alnreg_setSAM(opt, bns, pac, &s[i], p, 0x800 | 0x40<<i | extra_flag);
        mem_alnreg_formatSAM(opt, bns, &str, &s[i], p, NULL, 0);
        // g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p);
        // g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
        // aa[i][n_aa[i]++] = g[i];
      }
      s[i].sam = strdup(str.s); str.l = 0;
    }

    /* for (i = 0; i < n_aa[0]; ++i) */
    /*   mem_aln2sam(opt, bns, &str, &s[0], n_aa[0], aa[0], i, &h[1]); // write read1 hits  */
    /* s[0].sam = strdup(str.s); str.l = 0; */
    /*  */
    /* for (i = 0; i < n_aa[1]; ++i) */
    /*   mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]); // write read2 hits  */
    /* s[1].sam = str.s; */
    // if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);

    // Sanity check read names
    check_paired_read_names(s[0].name, s[1].name);

    // free CIGAR and XA
    for (i = 0; i < 2; ++i) {
      free(h[i].cigar); free(g[i].cigar);
      if (XA[i] == 0) continue;
      for (j = 0; j < a[i].n; ++j) free(XA[i][j]);
      free(XA[i]);
    }

  } else goto NO_PAIRING;
  return n;

NO_PAIRING:

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

