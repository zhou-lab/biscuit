/* utility functions for biscuit
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2023-2024 Jacob.Morrison@vai.org
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
**/

#include "bisc_utils.h"

const char nt256int8_to_methcode[NSTATUS_METH] = "RCN";     // R: retention, C: conversion, N: NA
const char nt256int8_to_basecode[NSTATUS_BASE] = "ACGTNYR"; // Y is [CT] and R is [GA]
const char *cytosine_context[]      = {"CG" ,"CHG" ,"CHH" ,"CG" ,"CHG","CHH","CN"};
const char *cytosine_context_nome[] = {"HCG","HCHG","HCHH","GCG","GCH","GCH","CN"};

cytosine_context_t fivenuc_context(refcache_t *rs, uint32_t rpos, char rb, char *fivenuc) {
    /* get five nucleotide context sequence */
    if (rpos == 1) { /* marginal cases, beginning of chromosome */
        subseq_refcache2(rs, 1, fivenuc+2, 3);
        fivenuc[0] = fivenuc[1] = 'N';
    } else if (rpos == 2) {
        subseq_refcache2(rs, 1, fivenuc+1, 4);
        fivenuc[0] = 'N';
    } else if (rpos == (unsigned) rs->seqlen) {  /* end of chromosome */
        subseq_refcache2(rs, rpos-2, fivenuc, 3);
        fivenuc[3] = fivenuc[4] = 'N';
    } else if (rpos == (unsigned) rs->seqlen-1) {
        subseq_refcache2(rs, rpos-2, fivenuc, 4);
        fivenuc[4] = 'N';
    } else {
        subseq_refcache2(rs, rpos-2, fivenuc, 5);
    }

    // reverse complement on 'G'
    if (rb == 'G') nt256char_rev_ip(fivenuc, 5);

    int i;
    for (i=0; i<5; ++i)
        if (fivenuc[i] == 'N')
            return CTXT_NA;

    if (rb != 'C' && rb != 'G') return CTXT_NA;

    /* classify into cytosine context classes */
    if (fivenuc[3] == 'G') {
        if (fivenuc[1]=='G') return CTXT_GCG;
        else return CTXT_HCG;
    } else if (fivenuc[4] == 'G') {
        if (fivenuc[1]=='G') return CTXT_GCHG;
        else return CTXT_HCHG;
    } else {
        if (fivenuc[1]=='G') return CTXT_GCHH;
        else return CTXT_HCHH;
    }
}

/* return -1 if abnormal (missing bsstrand) */
/* TODO: stratify by sequence context */
uint32_t cnt_retention(refcache_t *rs, bam1_t *b, uint8_t bsstrand) {
    uint32_t cnt = 0;

    bam1_core_t *c = &b->core;
    uint32_t i, rpos = c->pos+1, qpos = 0;
    uint32_t op, oplen;
    char rb, qb;
    unsigned j;
    for (i=0; i<c->n_cigar; ++i) {
        op = bam_cigar_op(bam_get_cigar(b)[i]);
        oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
        switch(op) {
            case BAM_CMATCH:
                for (j=0; j<oplen; ++j) {
                    rb = refcache_getbase_upcase(rs, rpos+j);
                    qb = bscall(b, qpos+j);
                    if (bsstrand) {
                        if (rb == 'C' && qb == 'C') cnt++;
                    } else {
                        if (rb == 'G' && qb == 'G') cnt++;
                    }
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case BAM_CINS:
                qpos += oplen;
                break;
            case BAM_CDEL:
                rpos += oplen;
                break;
            case BAM_CSOFT_CLIP:
                qpos += oplen;
                break;
            case BAM_CHARD_CLIP:
                qpos += oplen;
                break;
            default:
                fprintf(stderr, "Unknown cigar, %u\n", op);
                abort();
        }
    }

    return cnt;
}

uint32_t get_mate_length(char *m_cigar) {
    // reference length of mate read (excludes clipping and insertions)
    uint32_t length = 0;

    // An unmapped read will have a '*' as its CIGAR string
    // Treat these reads as having length 0 (this is why the initialization value is 0)
    if (*m_cigar != '*') {
        char *query;
        size_t n_cigar = 0;

        // Count number of CIGAR operations
        for (query = m_cigar; *m_cigar && *m_cigar != '\0'; ++m_cigar) {
            if (!isdigit((unsigned char) *m_cigar)) { ++n_cigar; }
        }

        if (*m_cigar++ != '\0') { wzfatal("Malformed MC tag CIGAR string\n"); }
        if (n_cigar == 0)       { wzfatal("No CIGAR operations found in MC tag\n"); }
        if (n_cigar >= 65536)   { wzfatal("Too many CIGAR operations found in MC tag\n"); }

        uint32_t i;
        int op;
        uint32_t *cigar = malloc(n_cigar * sizeof(uint32_t));
        for (i = 0; i < n_cigar; ++i, ++query) {
            cigar[i] = strtol(query, &query, 10) << BAM_CIGAR_SHIFT;

            op = (uint8_t)*query >= 128? -1 : bam_cigar_table[(int)*query];
            if (op < 0) { wzfatal("Unrecognized CIGAR operator\n"); }

            cigar[i] |= op;
        }

        length = bam_cigar2rlen(n_cigar, cigar);

        free(cigar);
    }

    return length;
}

uint8_t infer_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual) {

    /* infer bsstrand from nC2T and nG2A on high quality bases */

    bam1_core_t *c = &b->core;
    uint32_t i, rpos = c->pos+1, qpos = 0;
    uint32_t op, oplen;
    char rb, qb;
    int nC2T=0, nG2A=0; unsigned j;
    for (i=0; i<c->n_cigar; ++i) {
        op = bam_cigar_op(bam_get_cigar(b)[i]);
        oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
        switch(op) {
            case BAM_CMATCH:
                for (j=0; j<oplen; ++j) {
                    rb = refcache_getbase_upcase(rs, rpos+j);
                    qb = bscall(b, qpos+j);
                    if (bam_get_qual(b)[qpos+j] < min_base_qual) continue;
                    if (rb == 'C' && qb == 'T') nC2T++;
                    if (rb == 'G' && qb == 'A') nG2A++;
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case BAM_CINS:
                qpos += oplen;
                break;
            case BAM_CDEL:
                rpos += oplen;
                break;
            case BAM_CSOFT_CLIP:
                qpos += oplen;
                break;
            case BAM_CHARD_CLIP:
                qpos += oplen;
                break;
            default:
                fprintf(stderr, "Unknown cigar, %u\n", op);
                abort();
        }
    }
    if (nC2T >= nG2A) return 0;
    else return 1;
}

uint8_t get_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual, int allow_u) {
    uint8_t *s;

    /* bwa-meth flag has highest priority */
    s = bam_aux_get(b, "YD");
    if (s) {
        s++;
        if (*s == 'f') return 0;
        else if (*s == 'r') return 1;
        else if (*s == 'u' && allow_u) return 2;
    }

    /* bsmap flag */
    s = bam_aux_get(b, "ZS");
    if (s) {
        s++;
        if (*s == '+') return 0;
        else if (*s == '-') return 1;
    }

    /* bismark flag */
    s = bam_aux_get(b, "XG");
    if (s) {
        s++;
        if (strcmp((char*)s, "CT")==0) return 0;
        else if (strcmp((char*)s, "GA")==0) return 1;
    }

    /* otherwise, guess the bsstrand from nCT and nGA */
    return infer_bsstrand(rs, b, min_base_qual);
}

void pop_record_by_block_id(record_v *records, int64_t block_id, record_t *record) {
    uint64_t i;
    record_t *r;
    for (i=0; i<records->size; ++i) {
        r = ref_record_v(records, i);
        if (r->block_id == block_id) {
            *record = *r;             /* copy the record and set slot on shelf to OBSOLETE */
            r->block_id = RECORD_SLOT_OBSOLETE;
            return;
        }
    }
    record->block_id = RECORD_SLOT_OBSOLETE;
}

void put_into_record_v(record_v *records, record_t rec) {
    uint64_t i;
    record_t *r;

    /* fill blanks */
    for (i=0; i<records->size; ++i) {
        r = ref_record_v(records, i);
        if (r->block_id == RECORD_SLOT_OBSOLETE) {
            *r = rec;
            return;
        }
    }

    /* get a new slot */
    r = next_ref_record_v(records);
    *r = rec;
    return;
}

