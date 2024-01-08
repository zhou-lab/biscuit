/* multi-way pileup and somatic calling on cytosine converted bams
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021-2024 Jacob.Morrison@vai.org
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

#include <errno.h>
#include "pileup.h"

typedef struct {
    int n; /* number of sites */
    pileup_data_v **data;
} pileup_t;

pileup_t *init_pileup(int n) {
    pileup_t *p = malloc(sizeof(pileup_t));
    p->n = n;
    p->data = calloc(n, sizeof(pileup_data_v*));
    return p;
}

void destroy_pileup(pileup_t *p) {
    int i;
    for (i=0; i<p->n; ++i) {
        if (p->data[i]) free_pileup_data_v(p->data[i]);
    }
    free(p->data);
    free(p);
}

typedef struct {
    char **bam_fns;
    int n_bams;
    char *ref_fn;
    wqueue_t(window) *q;
    wqueue_t(record) *rq;
    pileup_conf_t *conf;
} result_t;

typedef struct {
    wqueue_t(record) *q;
    int n_bams;
    char **bam_fns;
    char *outfn;
    char *statsfn;
    char *header;
    target_v *targets;
    pileup_conf_t *conf;
} writer_conf_t;

static void print_meth_average_1chrom(FILE *out, char *sample, char *chrom, double *betasum, int64_t *cnt, writer_conf_t *c) {
    if (c->conf->comm.is_nome) { // nome-seq

        int64_t k_hcg  = cnt[CTXT_HCG];
        double  b_hcg  = betasum[CTXT_HCG];
        int64_t k_hchg = cnt[CTXT_HCHG];
        double  b_hchg = betasum[CTXT_HCHG];
        int64_t k_hchh = cnt[CTXT_HCHH];
        double  b_hchh = betasum[CTXT_HCHH];
        int64_t k_hch  = k_hchg + k_hchh;
        double  b_hch  = b_hchg + b_hchh;
        int64_t k_gch  = cnt[CTXT_GCG] +
            cnt[CTXT_GCHG] +
            cnt[CTXT_GCHH];
        double  b_gch  = betasum[CTXT_GCG] + 
            betasum[CTXT_GCHG] +
            betasum[CTXT_GCHH];

        if (k_hcg > 0) {             // skip chrom with no base coverage
            fprintf(out, "%s\t%s", sample, chrom);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_hcg,  b_hcg  / (double) k_hcg  * 100);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_hchg, b_hchg / (double) k_hchg * 100);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_hchh, b_hchh / (double) k_hchh * 100);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_hch,  b_hch  / (double) k_hch  * 100);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_gch,  b_gch  / (double) k_gch  * 100);
            fputc('\n', out);
        }

    } else {

        int64_t k_cg  = cnt[CTXT_GCG]+cnt[CTXT_HCG];
        double  b_cg  = betasum[CTXT_GCG]+betasum[CTXT_HCG];
        int64_t k_chg = cnt[CTXT_GCHG]+cnt[CTXT_HCHG];
        double  b_chg = betasum[CTXT_GCHG]+betasum[CTXT_HCHG];
        int64_t k_chh = cnt[CTXT_GCHH]+cnt[CTXT_HCHH];
        double  b_chh = betasum[CTXT_GCHH]+betasum[CTXT_HCHH];
        int64_t k_ch  = k_chg + k_chh;
        double  b_ch  = b_chg + b_chh;

        if (k_cg > 0) {             // skip chrom with no base coverage
            fprintf(out, "%s\t%s", sample, chrom);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_cg,  b_cg  / (double) k_cg  * 100);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_chg, b_chg / (double) k_chg * 100);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_chh, b_chh / (double) k_chh * 100);
            fprintf(out, "\t%"PRId64"\t%1.3f%%", k_ch,  b_ch  / (double) k_ch  * 100);
            fputc('\n', out);
        }
    }
}

static void print_meth_average1(FILE *out, char *sample, double *betasum, int64_t *cnt, writer_conf_t *c) {

    // genome-wide
    double  betasum0[NCONTXTS] = {0.0};
    int64_t cnt0[NCONTXTS] = {0};
    int k0;

    // by chromosome meth average
    unsigned k, t;
    for (k=0; k<c->targets->size; ++k) {
        t = get_target_v(c->targets, k).tid;
        print_meth_average_1chrom(out, sample, get_target_v(c->targets,t).name,
                betasum+k*NCONTXTS, cnt+k*NCONTXTS, c);
        // genome-wide
        for (k0=0; k0<NCONTXTS; ++k0) {
            cnt0[k0] += cnt[k*NCONTXTS+k0];
            betasum0[k0] += betasum[k*NCONTXTS+k0];
        }
    }

    // whole genome meth average
    print_meth_average_1chrom(out, sample, "WholeGenome", betasum0, cnt0, c);
}

void *write_func(void *data) {
    writer_conf_t *c = (writer_conf_t*) data;

    FILE *out;
    if (c->outfn) {
        out=fopen(c->outfn, "w");
        if (!out) {
            fprintf(stderr, "[%s:%d] Cannot open output file: %s\nAbort.\n", __func__, __LINE__, c->outfn);
            fflush(stderr);
            exit(1);
        }
    } else out=stdout;

    if (c->header) fputs(c->header, out);
    int64_t next_block = 0;
    record_v *records = init_record_v(20);

    /* statistics */
    int smpl_block = c->targets->size * NCONTXTS;
    int chrm_block = NCONTXTS;
    double *betasum=(double*)calloc(c->n_bams*smpl_block, sizeof(double));
    int64_t *cnt=(int64_t*)calloc(c->n_bams*smpl_block, sizeof(int64_t));
    while (1) {
        record_t rec;
        wqueue_get(record, c->q, &rec);
        if(rec.block_id == RECORD_QUEUE_END) break;
        if (rec.block_id == next_block) {
            do {
                if (rec.s.s) {
                    if (fputs(rec.s.s, out) < 0 && errno == EPIPE) exit(1);

                    /* statistics */
                    int sid;
                    for (sid=0; sid<c->n_bams; ++sid) {
                        /* meth average */
                        int i;
                        for (i=0; i<NCONTXTS; ++i) {
                            betasum[sid*smpl_block+rec.tid*chrm_block+i] += rec.betasum_context[sid*NCONTXTS+i];
                            cnt[sid*smpl_block+rec.tid*chrm_block+i] += rec.cnt_context[sid*NCONTXTS+i];
                        }
                    }
                }
                free(rec.s.s);
                free(rec.betasum_context);
                free(rec.cnt_context);

                /* get next block from shelf if available else return OBSOLETE 
                   and retrieve new block from queue  */
                next_block++;
                pop_record_by_block_id(records, next_block, &rec);
            } while (rec.block_id != RECORD_SLOT_OBSOLETE);
        } else {                    /* shelf the block if not next */
            put_into_record_v(records, rec);
        }
    }

    if ((!c->statsfn) && c->outfn) c->statsfn = strdup(c->outfn);

    if (c->statsfn) {

        // print meth level average
        char *outfn = calloc(strlen(c->statsfn)+20,1);
        strcpy(outfn, c->statsfn);
        strcat(outfn, "_meth_average.tsv");
        FILE *out = fopen(outfn, "w");

        if (c->conf->comm.is_nome)
            fprintf(out, "sample\tchrm\tHCGn\tHCGb\tHCHGn\tHCHGb\tHCHHn\tHCHHb\tHCHn\tHCHb\tGCn\tGCb\n");
        else
            fprintf(out, "sample\tchrm\tCGn\tCGb\tCHGn\tCHGb\tCHHn\tCHHb\tCHn\tCHb\n");

        int sid;
        for (sid=0; sid<c->n_bams; ++sid)
            print_meth_average1(out, c->bam_fns[sid], betasum+sid*smpl_block, cnt+sid*smpl_block, c);

        fclose(out);
        free(outfn);
    }

    free(c->statsfn);
    free(cnt);
    free(betasum);

    free_record_v(records);
    if (c->outfn) {    /* for stdout, will close at the end of main */
        fflush(out);
        fclose(out);
    }
    return 0;
}

static void verbose_format(uint8_t bsstrand, pileup_data_v *dv, kstring_t *s, int sid) {

    uint32_t i, nf;

    /* return if no record match the bsstrand */
    int n=0;
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv,i);
        if (d->sid != sid) continue;
        if (d->bsstrand == bsstrand) ++n;
    }
    if (!n) return;

    char b='0'+bsstrand;

    /* 1. base */
    ksprintf(s, ";Bs%c=", b);
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv,i);
        if (d->sid != sid) continue;
        if (d->bsstrand == bsstrand) kputc(d->qb, s);
    }

    /* 2. status array */
    ksprintf(s, ";Sta%c=", b);
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv,i);
        if (d->sid != sid) continue;
        if (d->bsstrand == bsstrand) kputc('0'+(d->stat&0xf), s);
    }

    /* 3. base quality */
    ksprintf(s, ";Bq%c=", b);
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv,i);
        if (d->sid != sid) continue;
        if (d->bsstrand == bsstrand) kputc(d->qual+33, s);
    }

    /* 4. strand */
    ksprintf(s, ";Str%c=", b);
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv,i);
        if (d->sid != sid) continue;
        if (d->bsstrand == bsstrand) kputc(d->strand?'-':'+', s);
    }

    /* 5. position on read */
    ksprintf(s, ";Pos%c=", b);
    nf = 0;
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv,i);
        if (d->sid != sid) continue;
        if (d->bsstrand == bsstrand) {
            if (nf) kputc(',', s);
            else nf = 1;
            kputuw(d->qpos, s);
        }
    }

    /* 6. retention count 
       retention count, for diagnosing incomplete converted
       reads from CpH sites and mitochondrial sites */
    ksprintf(s, ";Rret%c=", b);
    nf = 0;
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv,i);
        if (d->sid != sid) continue;
        if (d->bsstrand == bsstrand) {
            if (nf) kputc(',', s);
            else nf = 1;
            kputuw(d->cnt_ret, s);
        }
    }
}

static int top_mutant(int *cnts_base1, uint32_t supp[NSTATUS_BASE], uint8_t rb_code) {

    uint8_t i;
    for (i=0; i<NSTATUS_BASE; ++i) {
        if (i!=BASE_N) supp[i] = (cnts_base1[i]<<4) | i;
        else supp[i] = 0;
    }

    qsort(supp, NSTATUS_BASE, sizeof(uint32_t), compare_supp);

    int cm1 = -1; // top mutant for genotyping
    for (i=0; i<NSTATUS_BASE; ++i) {
        uint8_t base = supp[i]&0xf;
        if (base == BASE_R && (rb_code == BASE_A || rb_code == BASE_G)) continue;
        if (base == BASE_Y && (rb_code == BASE_C || rb_code == BASE_T)) continue;
        if (base != BASE_N && base!=rb_code && (supp[i]>>4)>0) {
            cm1 = base;
            break;
        }
    }
    return cm1;
}

// redistribute base supports, considering bisulfite treatment
// special treatment here: (should have an option to turn this off)
// if (reference shows T or there is evidence of T) and there is no sign of C,
// we assume all counts of Y is T, vice versa and similar applies to G/A and R
static void redistribute_cnts(int *cnts_base, int n_bams, uint8_t rb_code) {

    int cnts_base_all[NSTATUS_BASE] = {0}; // all samples, this allows cross-sample inference
    int sid; uint8_t i;
    for (sid=0; sid<n_bams; ++sid) {
        for (i=0; i<NSTATUS_BASE; ++i)
            cnts_base_all[i] += cnts_base[sid*NSTATUS_BASE+i];
    }

    for (sid=0; sid<n_bams; ++sid) {
        int *cnts_base1 = cnts_base + NSTATUS_BASE*sid;
        if (((rb_code == BASE_T || cnts_base_all[BASE_T]) && cnts_base_all[BASE_C] == 0 && rb_code != BASE_C)) {
            cnts_base1[BASE_T] += cnts_base1[BASE_Y];
            cnts_base1[BASE_Y] = 0;
        }

        if (((rb_code == BASE_C || cnts_base_all[BASE_C]) && cnts_base_all[BASE_T] == 0 && rb_code != BASE_T)) {
            cnts_base1[BASE_C] += cnts_base1[BASE_Y];
            cnts_base1[BASE_Y] = 0;
        }

        if (((rb_code == BASE_A || cnts_base_all[BASE_A]) && cnts_base_all[BASE_G] == 0 && rb_code != BASE_G)) {
            cnts_base1[BASE_A] += cnts_base1[BASE_R];
            cnts_base1[BASE_R] = 0;
        }

        if (((rb_code == BASE_G || cnts_base_all[BASE_G]) && cnts_base_all[BASE_A] == 0 && rb_code != BASE_A)) {
            cnts_base1[BASE_G] += cnts_base1[BASE_R];
            cnts_base1[BASE_R] = 0;
        }
    }
}
  
static void plp_getcnts(pileup_data_v *dv, pileup_conf_t *conf, int *cnts_meth, int *cnts_base) {

    if (!dv) return;

    uint32_t i;
    for (i=0; i<dv->size; ++i) {
        pileup_data_t *d = ref_pileup_data_v(dv, i);
        /* read-position-based filtering */
        if (d->qual < conf->filt.min_base_qual) continue;
        // TODO: On the 3'-end, it should be before the end of mapping instead of
        // end of read since adaptor sequences are soft-clipped.
        if (d->qpos <= conf->filt.min_dist_end_5p || d->rlen < d->qpos + conf->filt.min_dist_end_3p) continue;
        cnts_meth[d->sid * NSTATUS_METH + (d->stat&0xf)]++;
        cnts_base[d->sid * NSTATUS_BASE + (d->stat>>4)]++;
    }
}

void pileup_genotype(int cref, int altsupp, pileup_conf_t *conf, char gt[4], double *_gl0, double *_gl1, double *_gl2, double *_gq) {

    double gl0=-1, gl1=-1, gl2=-1, gq=-1;
    if (cref >=0 || altsupp >= 0) {
        gl0 = log(conf->prior0) + genotype_lnlik(HOMOREF, cref, altsupp, conf->error, conf->contam);
        gl1 = log(conf->prior1) + genotype_lnlik(HET, cref, altsupp, conf->error, conf->contam);
        gl2 = log(conf->prior2) + genotype_lnlik(HOMOVAR, cref, altsupp, conf->error, conf->contam);
        if (gl0>gl1) {
            if (gl0>gl2) {
                gq = pval2qual(1 - exp(gl0 - ln_sum3(gl0, gl1, gl2)));
                strcpy(gt, "0/0");
            } else {
                gq = pval2qual(1 - exp(gl2 - ln_sum3(gl0, gl1, gl2)));
                strcpy(gt, "1/1");
            }
        } else if (gl1>gl2) {
            gq = pval2qual(1 - exp(gl1 - ln_sum3(gl0, gl1, gl2)));
            strcpy(gt, "0/1");
        } else {
            gq = pval2qual(1 - exp(gl2 - ln_sum3(gl0, gl1, gl2)));
            strcpy(gt, "1/1");
        }
    }
    *_gl0 = gl0; *_gl1 = gl1; *_gl2 = gl2; *_gq = gq;
}

static void plp_format(refcache_t *rs, char *chrm, uint32_t rpos, pileup_data_v *dv, pileup_conf_t *conf,
        record_t *rec, int n_bams, int *cnts_meth, int *cnts_base, int *cnts_base_redist) {

    kstring_t *s = &rec->s;

    uint32_t i;
    char rb = refcache_getbase_upcase(rs, rpos);
    if (rb == 'N') return;
    int rb_code = nt256char_to_nt256int8_table[(uint8_t)rb];

    memset(cnts_meth, 0, sizeof(int)*NSTATUS_METH*n_bams);
    memset(cnts_base, 0, sizeof(int)*NSTATUS_BASE*n_bams);
    plp_getcnts(dv, conf, cnts_meth, cnts_base);

    memcpy(cnts_base_redist, cnts_base, sizeof(int)*NSTATUS_BASE*n_bams);
    if (conf->ambi_redist) redistribute_cnts(cnts_base_redist, n_bams, rb_code);

    int cnts_base_allsamples[NSTATUS_BASE] = {0};
    int cnts_meth_allsamples[NSTATUS_METH] = {0};
    int sid;
    for (sid=0; sid<n_bams; ++sid) {
        for (i=0; i<NSTATUS_METH; ++i) cnts_meth_allsamples[i] += cnts_meth[sid*NSTATUS_METH+i];
        for (i=0; i<NSTATUS_BASE; ++i) cnts_base_allsamples[i] += cnts_base_redist[sid*NSTATUS_BASE+i];
    }

    uint32_t supp[NSTATUS_BASE];
    int cm1 = top_mutant(cnts_base_allsamples, supp, rb_code);

    /* if not SNP but no signal for METH_RETENTION or METH_CONVERSION,
       skip the print when in non-verbose mode */
    if (cm1 < 0 && !conf->comm.verbose
            && cnts_meth_allsamples[METH_RETENTION] == 0
            && cnts_meth_allsamples[METH_CONVERSION] == 0)
        return;

    /* initialize genotyping */
    char **gt = malloc(n_bams*sizeof(char*));
    for (sid=0; sid<n_bams; ++sid) {
        gt[sid] = malloc(4*sizeof(char));
        strcpy(gt[sid], "./.");
    }
    double *gl0 = malloc(n_bams*sizeof(double));
    double *gl1 = malloc(n_bams*sizeof(double));
    double *gl2 = malloc(n_bams*sizeof(double));
    double *gq = malloc(n_bams*sizeof(double));
    /* genotype each sample */
    int any_methcallable = 0;
    double lowest_gq = 0;
    uint8_t *methcallable = calloc(n_bams, sizeof(uint8_t));
    for (sid=0; sid<n_bams; ++sid) {
        int *cnts_base1 = cnts_base_redist + NSTATUS_BASE*sid;
        int *cnts_meth1 = cnts_meth + NSTATUS_METH*sid;

        /* determine if SNP is interfering methylation calling */
        if (cnts_meth1[METH_RETENTION] + cnts_meth1[METH_CONVERSION] > 0) {
            //if (cnts_base1[BASE_T]==0 && rb == 'C') methcallable[sid] = 1;
            //if (cnts_base1[BASE_A]==0 && rb == 'G') methcallable[sid] = 1;
            if (rb == 'C') {
                if (cnts_base1[BASE_T] == 0)
                    methcallable[sid] = 1;
                else if (cnts_base1[BASE_C]>0 && cnts_base1[BASE_T]/(double) cnts_base1[BASE_C] < 0.05)
                    methcallable[sid] = 1;
            }

            if (rb == 'G') {
                if (cnts_base1[BASE_A] == 0)
                    methcallable[sid] = 1;
                else if (cnts_base1[BASE_G]>0 && cnts_base1[BASE_A]/(double) cnts_base1[BASE_G] < 0.05)
                    methcallable[sid] = 1;
            }
        }

        gl0[sid] = -1; gl1[sid] = -1; gl2[sid] = -1; gq[sid] = 0;
        int nref = cnts_base1[rb_code];
        int nalt = cm1 >= 0 ? cnts_base1[cm1] : 0;
        if (nref + nalt > 0) pileup_genotype(nref, nalt, conf, gt[sid], gl0+sid, gl1+sid, gl2+sid, gq+sid);

        if (gq[sid] < lowest_gq || !sid) lowest_gq = gq[sid];

        if (methcallable[sid])
            any_methcallable = 1;
    }

    /* determine somatic mutation status */
    double squal=0.0; int ss=5;
    if (conf->somatic && cm1>=0) {
        uint32_t supp_t[NSTATUS_BASE];
        int cm1_t = top_mutant(cnts_base_redist, supp_t, rb_code);

        if (cm1_t >= 0) {
            int altcnt_t = cnts_base_redist[cm1_t];
            int altcnt_n = cnts_base_redist[NSTATUS_BASE+cm1_t];
            int cref_t = cnts_base_redist[rb_code];
            int cref_n = cnts_base_redist[NSTATUS_BASE+rb_code];
            squal = pval2qual(somatic_posterior(cref_t, altcnt_t, cref_n, altcnt_n, conf->error, conf->mu, conf->mu_somatic, conf->contam));
            if (squal>1) ss=2;
            else if (gt[1][2]=='1') {
                ss=1;
            } else {
                ss=0;
            }
        }
    }

    /****** output *****/
    /* CHROM, POS, ID, REF */
    ksprintf(s, "%s\t%u\t.\t%c\t", chrm, rpos, rb);

    /* ALT */
    if (cm1 >= 0) {
        char m;
        if (cm1 == BASE_Y || cm1 == BASE_R) m = 'N';
        else m = nt256int8_to_basecode[cm1];
        kputc(m, s);
    } else {
        kputc('.', s);
    }

    /* QUAL */
    ksprintf(s, "\t%d", (int) lowest_gq);
    if (lowest_gq > 5) {
        kputs("\tPASS\t", s);
    } else {
        kputs("\tLowQual\t", s);
    }

    /* INFO tags */
    cytosine_context_t ctt=CTXT_NA;
    ksprintf(s, "NS=%d", n_bams);
    if (rb == 'C' || rb == 'G') {
        char fivenuc[5];
        ctt = fivenuc_context(rs, rpos, rb, fivenuc);
        ksprintf(s, ";CX=%s", conf->comm.is_nome?cytosine_context_nome[ctt]:cytosine_context[ctt]);
        ksprintf(s, ";N5=%.5s", fivenuc);
    }
    if (conf->somatic && cm1>=0) {
        ksprintf(s, ";SS=%d", ss);
        ksprintf(s, ";SC=%d", (int) squal);
    }
    if (cm1 >= 0 && (cm1 == BASE_Y || cm1 == BASE_R)) {
        kputs(";AB=", s);
        kputc(nt256int8_to_basecode[cm1], s);
    }

    /* FORMAT */
    kputs("\tGT:GL1:GQ:DP", s);
    kputs(":SP", s);
    if (cm1 >= 0) kputs(":AC:AF1", s);
    if (any_methcallable) kputs(":CV:BT", s);

    /* loop over samples print format */
    double beta;
    for (sid=0; sid<n_bams; ++sid) {
        int *cnts_base1 = cnts_base + NSTATUS_BASE*sid;
        int *cnts_base1_redist = cnts_base_redist + NSTATUS_BASE*sid;
        int *cnts_meth1 = cnts_meth + NSTATUS_METH*sid;
        int dp=0;
        if (dv) for(i=0; i<dv->size; ++i)	if (ref_pileup_data_v(dv, i)->sid == sid) ++dp;

        /* GT, GL, GQ */
        if (gq[sid]>0 && dp) {
            ksprintf(s, "\t%s:%1.0f,%1.0f,%1.0f:%1.0f", gt[sid], 
                    max(-1000, gl0[sid]), max(-1000, gl1[sid]), max(-1000, gl2[sid]), gq[sid]);
        } else {
            ksprintf(s, "\t./.:.,.,.:0");
        }

        /* DP */
        if (dp) ksprintf(s, ":%d", dp);
        else kputs(":0", s);

        /* SP */
        kputc(':', s);
        int added = 0;
        if (cnts_base1[rb_code]) { ksprintf(s, "%c%d", rb, cnts_base1[rb_code]); added = 1; }
        uint8_t i;
        for (i=0; i<NSTATUS_BASE; ++i) {
            if (i == BASE_N) continue;
            if (i == rb_code) continue;
            if (cnts_base1[i] <= 0) continue;
            ksprintf(s, "%c%d", nt256int8_to_basecode[i], cnts_base1[i]);
            added = 1;
        }
        if (!added) kputc('.', s);

        /* AC, AF */
        if (cm1 >= 0) {
            int nref = cnts_base1_redist[rb_code];
            int nalt = cnts_base1_redist[cm1];
            ksprintf(s, ":%d:", nref+nalt);
            if (nref+nalt) ksprintf(s, "%1.2f", nalt/(double) (nref+nalt));
            else kputc('.', s);
        }

        /* CV, BT */
        if (any_methcallable) {
            if (methcallable[sid]) {
                beta = (double) cnts_meth1[METH_RETENTION] / (double) (cnts_meth1[METH_RETENTION]+cnts_meth1[METH_CONVERSION]);
                if (ctt != CTXT_NA) {
                    rec->betasum_context[sid*NCONTXTS+ctt] += beta;
                    rec->cnt_context[sid*NCONTXTS+ctt]++;
                }
                ksprintf(s, ":%d:%1.3f", cnts_meth1[METH_RETENTION]+cnts_meth1[METH_CONVERSION], beta);
            } else {
                kputs(":0:.", s);
            }
        }

        /* additional information printed on verbose theoretically this should be 
           put to FORMAT since they are sample-specific, but in FORMAT is usually 
           harder to parse, so they are appended to INFO these are not intended for
           formal submission, just for diagnostic purposes.  */
        if (conf->comm.verbose) {
            kputs("\tDIAGNOSE", s);
            if (methcallable)
                ksprintf(s, ";RN=%d;CN=%d", cnts_meth1[METH_RETENTION], cnts_meth1[METH_CONVERSION]);
            verbose_format(0, dv, s, sid);
            verbose_format(1, dv, s, sid);
        }
    }

    kputc('\n', s);
    for (sid=0; sid<n_bams; ++sid) free(gt[sid]);
    free(gt); free(gl0); free(gl1); free(gl2); free(gq);
    free(methcallable);
}

static void *process_func(void *_result) {

    result_t *res = (result_t*) _result;
    pileup_conf_t *conf = (pileup_conf_t*) res->conf;

    /* TODO: protect against res->n_bams == 0 */

    // open bam files
    int sid=0;                    // sample id
    htsFile **in_fhs = calloc(res->n_bams, sizeof(htsFile*));
    hts_idx_t **idxs = calloc(res->n_bams, sizeof(hts_idx_t*));
    for (sid=0; sid<res->n_bams; ++sid) {
        in_fhs[sid] = hts_open(res->bam_fns[sid], "rb");
        if (!in_fhs[sid]) {
            fprintf(stderr, "[%s:%d] Cannot open %s\nAbort.\n", __func__, __LINE__, res->bam_fns[sid]);
            fflush(stderr);
            exit(1);
        }
        idxs[sid] = sam_index_load(in_fhs[sid], res->bam_fns[sid]);
        if (!idxs[sid]) {
            fprintf(stderr, "[%s:%d] Cannot find index for %s\n", __func__, __LINE__, res->bam_fns[sid]);
            fflush(stderr);
            exit(1);
        }
    }
    /* header of first bam */
    // TODO: This assumes that the headers are the same for all input BAMs
    //       This is probably okay to assume, but it's better to use the respective headers for each BAM
    bam_hdr_t *bam_hdr = sam_hdr_read(in_fhs[0]);

    refcache_t *rs = init_refcache(res->ref_fn, 1000, 1000);
    uint32_t i; unsigned j;

    record_t rec;
    window_t w;
    while (1) {

        wqueue_get(window, res->q, &w);
        if (w.tid == -1) break;

        /* prepare record string */
        rec.s.l = rec.s.m = 0; rec.s.s = 0;
        rec.block_id = w.block_id;

        /* prepare statistics estimates */
        rec.tid = w.tid;
        rec.betasum_context = calloc(NCONTXTS*res->n_bams, sizeof(double));
        rec.cnt_context = calloc(NCONTXTS*res->n_bams, sizeof(int64_t));

        pileup_t *plp = init_pileup(w.end - w.beg);

        // chrm based on the first bam
        char *chrm = bam_hdr->target_name[w.tid];

        /* prepare reference */
        refcache_fetch(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);

        /* loop over bams */
        for (sid=0; sid<res->n_bams; ++sid) {
            htsFile *in = in_fhs[sid];
            hts_idx_t *idx = idxs[sid];

            hts_itr_t *iter = sam_itr_queryi(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
            bam1_t *b = bam_init1(); int ret;
            // loop over reads
            while ((ret = sam_itr_next(in, iter, b))>0) {

                uint8_t bsstrand = get_bsstrand(rs, b, conf->filt.min_base_qual, 0);

                /* read-based filtering */
                bam1_core_t *c = &b->core;
                if (c->qual < conf->filt.min_mapq) continue;
                if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->filt.min_read_len) continue;
                if (c->flag > 0){
                    if (conf->filt.filter_secondary && (c->flag & BAM_FSECONDARY)) continue;
                    if (conf->filt.filter_duplicate && (c->flag & BAM_FDUP)) continue;
                    if (conf->filt.filter_ppair && c->flag & BAM_FPAIRED && !(c->flag & BAM_FPROPER_PAIR)) continue;
                    if (conf->filt.filter_qcfail && c->flag & BAM_FQCFAIL) continue;
                }

                uint8_t *nm = bam_aux_get(b, "NM");
                if (nm && bam_aux2i(nm) > conf->filt.max_nm) continue;

                uint8_t *as = bam_aux_get(b, "AS");
                if (as && bam_aux2i(as) < conf->filt.min_score) continue;

                uint32_t cnt_ret = cnt_retention(rs, b, bsstrand);
                if (cnt_ret > conf->filt.max_retention) continue;

                // read positions and lengths
                uint32_t rpos = c->pos+1, qpos = 0, rmpos = c->mpos + 1;
                uint32_t read_length = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));

                uint32_t mate_length;
                uint8_t *mc = bam_aux_get(b, "MC");
                if (mc) {
                    mc++;
                    mate_length = get_mate_length((char *)mc);
                } else {
                    // If MC tag is missing, then assume reads are the same length
                    mate_length = read_length;
                }

                // -1 accounts for the 1-based nature of the coordinates
                uint32_t rend  = rpos  + read_length - 1;
                uint32_t rmend = rmpos + mate_length - 1;

                for (i=0; i<c->n_cigar; ++i) {
                    uint32_t op = bam_cigar_op(bam_get_cigar(b)[i]);
                    uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                    char qb, rb;
                    switch(op) {
                        case BAM_CMATCH:
                            for (j=0; j<oplen; ++j) {

                                if (rpos+j<w.beg || rpos+j>=w.end) continue; /* include begin but not end */
                                rb = refcache_getbase_upcase(rs, rpos+j);
                                qb = bscall(b, qpos+j);

                                /* If reads 1 and 2 overlap, skip counting bases in read 2
                                 * Read lengths are relative to the reference as defined by the CIGAR string (if MC tag
                                 *     not given in read, then mate length is assumed to be the same as the current
                                 *     read)
                                 * Overlapping bases in both proper and improper pairs (if user requests these be
                                 *     included) will be ignored
                                 *  
                                 * The filtering removes bases from read 2 (usually the read on the complement strand)
                                 * that fall into the overlapped region.
                                 */
                                if ((conf->filt.filter_doublecnt) &&
                                    (c->flag & BAM_FREAD2) &&
                                    (rpos+j >= max(rpos, rmpos)) &&
                                    (rpos+j <= min(rend, rmend))) {
                                    continue;
                                }

                                pileup_data_v **plp_data_vec = plp->data+rpos+j-w.beg;
                                if (!*plp_data_vec) *plp_data_vec = init_pileup_data_v(2);
                                pileup_data_t *d = next_ref_pileup_data_v(*plp_data_vec);
                                d->sid = sid;
                                d->qual = bam_get_qual(b)[qpos+j];
                                d->cnt_ret = (unsigned) cnt_ret;
                                d->strand = (c->flag&BAM_FREVERSE)?1:0;
                                d->qpos = qpos+j+1;
                                d->rlen = c->l_qseq;
                                d->bsstrand = bsstrand;
                                d->qb = qb;

                                d->stat = 0;
                                if (bsstrand) {	/* BSC */
                                    if (rb == 'G') {
                                        if (qb == 'A') d->stat = METH_CONVERSION;
                                        else if (qb == 'G') d->stat = METH_RETENTION;
                                        else d->stat = METH_NA;
                                    } else d->stat = METH_NA;
                                    if (qb == 'A') d->stat |= BASE_R << 4;
                                    else d->stat |= (nt256char_to_nt256int8_table[(uint8_t)qb]<<4);
                                } else {		/* BSW */
                                    if (rb == 'C') {
                                        if (qb == 'T') d->stat = METH_CONVERSION;
                                        else if (qb == 'C') d->stat = METH_RETENTION;
                                        else d->stat = METH_NA;
                                    } else d->stat = METH_NA;
                                    if (qb == 'T') d->stat |= BASE_Y << 4;
                                    else d->stat |= (nt256char_to_nt256int8_table[(uint8_t)qb]<<4);
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
            }

            bam_destroy1(b);
            hts_itr_destroy(iter);
        }

        /* loop over cytosines and format */
        int *cnts_meth = calloc(NSTATUS_METH*res->n_bams, sizeof(int));
        int *cnts_base = calloc(NSTATUS_BASE*res->n_bams, sizeof(int));
        int *cnts_base_redist = calloc(NSTATUS_BASE*res->n_bams, sizeof(int));
        for (j=w.beg; j<w.end; ++j) {
            pileup_data_v *plp_data = plp->data[j-w.beg];
            if (plp_data) {
                plp_format(rs, chrm, j, plp_data, res->conf, &rec, res->n_bams, cnts_meth, cnts_base, cnts_base_redist);
            }
        }
        free(cnts_meth); free(cnts_base); free(cnts_base_redist);

        /* put output string to output queue */
        wqueue_put2(record, res->rq, rec);

        destroy_pileup(plp);
    }
    free_refcache(rs);
    for (sid=0; sid<res->n_bams; ++sid) {
        hts_close(in_fhs[sid]);
        hts_idx_destroy(idxs[sid]);
    }
    free(in_fhs);
    free(idxs);
    bam_hdr_destroy(bam_hdr);
    return 0;
}

static void head_append_verbose(char *pb, char b, kstring_t *s) {
    ksprintf(s, "##FORMAT=<ID=Bs%c,Number=1,Type=String,Description=\"base identity, %s\">\n", b, pb);
    ksprintf(s, "##FORMAT=<ID=Sta%c,Number=1,Type=String,Description=\"Status code, %s (0,1,2 for retention, conversion and NA)\">\n", b, pb);
    ksprintf(s, "##FORMAT=<ID=Bq%c,Number=1,Type=String,Description=\"base quality, %s\">\n", b, pb);
    ksprintf(s, "##FORMAT=<ID=Str%c,Number=1,Type=String;Description=\"strands, %s\">\n", b, pb);
    ksprintf(s, "##FORMAT=<ID=Pos%c,Number=1,Type=String;Description=\"position in read, %s\">\n", b, pb);
    ksprintf(s, "##FORMAT=<ID=Rret%c,Number=1,Type=String;Description=\"Number of retention in read, %s\">\n", b, pb);
}

char *print_vcf_header(char *reffn, target_v *targets, char **argv, int argc, pileup_conf_t *conf, char **in_fns, int n_fns) {

    kstring_t header; header.l = header.m = 0; header.s = 0;
    kputs("##fileformat=VCFv4.1\n", &header);
    ksprintf(&header, "##reference=%s\n", reffn);
    ksprintf(&header, "##source=biscuitV%s\n", PACKAGE_VERSION);

    unsigned j;
    for (j=0; j<targets->size; ++j) {
        target_t *t = ref_target_v(targets, j);
        ksprintf(&header, "##contig=<ID=%s,length=%d>\n", t->name, t->len);
    }
    kputs("##program=<cmd=biscuit", &header);
    int i;
    for (i=0; i<argc; ++i)
        ksprintf(&header, " %s", argv[i]);
    kputs(">\n", &header);
    kputs("##FILTER=<ID=PASS,Description=\"All filters passed\">\n", &header);
    kputs("##FILTER=<ID=LowQual,Description=\"Genotype quality smaller than 5\">\n", &header);
    kputs("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n", &header);
    if (conf->comm.is_nome)
        kputs("##INFO=<ID=CX,Number=1,Type=String,Description=\"Cytosine context (HCG, HCHG, HCHH, GCG, GCH)\">\n", &header);
    else
        kputs("##INFO=<ID=CX,Number=1,Type=String,Description=\"Cytosine context (CG, CHH or CHG)\">\n", &header);
    kputs("##INFO=<ID=N5,Number=1,Type=String,Description=\"5-nucleotide context, centered around target cytosine\">\n", &header);
    kputs("##INFO=<ID=AB,Number=A,Type=String,Description=\"When true alt-allele is ambiguous, ALT field will be N and true alt-allele is stored here, following IUPAC code convention. This option does not appear when ALT != N.\">\n", &header);

    if (conf->somatic) {
        kputs("##INFO=<ID=SS,Number=1,Type=String,Description=\"Somatic status 0) WILDTYPE; 1) GERMLINE; 2) SOMATIC; 3) LOH; 4) POST_TRX_MOD; 5) UNKNOWN;\">\n", &header);
        kputs("##INFO=<ID=SC,Number=1,Type=Float,Description=\"Somatic score\">\n", &header);
        kputs("##INFO=<ID=AF1,Number=1,Type=Float,Description=\"Variant allele fraction\">\n", &header);
    }

    kputs("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n", &header);
    kputs("##FORMAT=<ID=SP,Number=.,Type=String,Description=\"Allele support (considering bisulfite conversion, with filtering)\">\n", &header);
    kputs("##FORMAT=<ID=AC,Number=.,Type=Integer,Description=\"Depth in calculating alternative allele frequency (after inference, with filtering)\">\n", &header);
    kputs("##FORMAT=<ID=AF1,Number=.,Type=Float,Description=\"Alternative allele frequency (after inference, with filtering)\">\n", &header);
    kputs("##FORMAT=<ID=CV,Number=1,Type=Integer,Description=\"Effective (strand-specific) coverage on cytosine\">\n", &header);
    kputs("##FORMAT=<ID=BT,Number=1,Type=Float,Description=\"Cytosine methylation fraction (aka beta value, with filtering)\">\n", &header);
    kputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype from normal\">\n", &header);
    kputs("##FORMAT=<ID=GL1,Number=3,Type=Float,Description=\"Genotype likelihoods for the first alternative allele\">\n", &header);
    kputs("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality (phred-scaled)\">\n", &header);

    if (conf->comm.verbose) {
        kputs("##FORMAT=<ID=RN,Number=1,Type=Integer,Description=\"Retention count (with filtering)\">\n", &header);
        kputs("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Conversion count (with filtering)\">\n", &header);
        char plpbsstrand[4];
        strcpy(plpbsstrand, "BSW");
        head_append_verbose(plpbsstrand, '0', &header);
        strcpy(plpbsstrand, "BSC");
        head_append_verbose(plpbsstrand, '1', &header);
    }
    kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", &header);

    // inference of sample name from bam file name
    int sid=0;
    for (sid=0; sid<n_fns; ++sid) {
        kputc('\t', &header);
        char *path=strdup(in_fns[sid]);
        char *bname = basename(path);
        if (strcmp(bname+strlen(bname)-4,".bam")==0)
            bname[strlen(bname)-4]=0;
        kputs(bname, &header);
        free(path);
    }
    kputc('\n', &header);

    return header.s;
}

void pileup_conf_init(pileup_conf_t *conf) {
    conf->comm = bisc_common_init();
    conf->bt = bisc_threads_init();
    conf->filt = meth_filter_init();

    conf->ambi_redist = 1;

    /* genotyping */
    conf->somatic = 0;
    conf->error = 0.001;
    conf->mu = 0.001;
    conf->mu_somatic = 0.001;
    conf->contam = 0.01;
    conf->prior1 = 0.33333;
    conf->prior2 = 0.33333;
    conf->prior0 = 1.0 - conf->prior1 - conf->prior2;
    if (conf->prior0 < 0) { fprintf(stderr, "[Error] genotype prior0 (%1.3f) must be from 0 to 1. \n", conf->prior0); exit(1); }
    if (conf->prior1 < 0) { fprintf(stderr, "[Error] genotype prior1 (%1.3f) must be from 0 to 1. \n", conf->prior1); exit(1); }
    if (conf->prior2 < 0) { fprintf(stderr, "[Error] genotype prior2 (%1.3f) must be from 0 to 1. \n", conf->prior2); exit(1); }
}

static int usage(pileup_conf_t *conf) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit pileup [options] <ref.fa> <in1.bam> [in2.bam in3.bam ...]\n");
    fprintf(stderr, "Som. Mode Usage: biscuit pileup [options] <-S -T tum.bam -I norm.bam> <ref.fa>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -g STR      Region (optional, will process the whole bam if not specified)\n");
    fprintf(stderr, "    -@ INT      Number of threads [%d]\n", conf->bt.n_threads);
    fprintf(stderr, "    -s INT      Step of window dispatching [%d]\n", conf->bt.step);
    fprintf(stderr, "    -N          NOMe-seq mode [off]\n");
    fprintf(stderr, "    -S          Somatic mode, must provide -T and -I arguments [off]\n");
    fprintf(stderr, "    -T STR      Somatic mode, tumor BAM\n");
    fprintf(stderr, "    -I STR      Somatic mode, normal BAM\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o STR      Output file [stdout]\n");
    fprintf(stderr, "    -w STR      Pileup statistics output prefix [same as output]\n");
    fprintf(stderr, "    -v INT      Verbosity level (0: no added info printed, 0<INT<=5: print\n");
    fprintf(stderr, "                    diagnostic info, INT>5: print diagnostic and debug info) [0]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -b INT      Minimum base quality [%u]\n", conf->filt.min_base_qual);
    fprintf(stderr, "    -m INT      Minimum mapping quality [%u]\n", conf->filt.min_mapq);
    fprintf(stderr, "    -a INT      Minimum alignment score (from AS-tag) [%u]\n", conf->filt.min_score);
    fprintf(stderr, "    -t INT      Maximum cytosine retention in a read [%u]\n", conf->filt.max_retention);
    fprintf(stderr, "    -l INT      Minimum read length [%u]\n", conf->filt.min_read_len);
    fprintf(stderr, "    -5 INT      Minimum distance to 5' end of a read [%u]\n", conf->filt.min_dist_end_5p);
    fprintf(stderr, "    -3 INT      Minimum distance to 3' end of a read [%u]\n", conf->filt.min_dist_end_3p);
    fprintf(stderr, "    -r          NO redistribution of ambiguous (Y/R) calls in SNP genotyping\n");
    fprintf(stderr, "    -c          NO filtering secondary mapping\n");
    fprintf(stderr, "    -d          Double count cytosines in overlapping mate reads (avoided\n");
    fprintf(stderr, "                    by default)\n");
    fprintf(stderr, "    -u          NO filtering of duplicate flagged reads\n");
    fprintf(stderr, "    -p          NO filtering of improper pair flagged reads\n");
    fprintf(stderr, "    -n INT      Maximum NM tag [%d]\n", conf->filt.max_nm);
    fprintf(stderr, "\n");
    fprintf(stderr, "Genotyping options:\n");
    fprintf(stderr, "    -E FLOAT    Error rate [%1.3f]\n", conf->error);
    fprintf(stderr, "    -M FLOAT    Mutation rate [%1.3f]\n", conf->mu);
    fprintf(stderr, "    -x FLOAT    Somatic mutation rate [%1.3f]\n", conf->mu_somatic);
    fprintf(stderr, "    -C FLOAT    Contamination rate [%1.3f]\n", conf->contam);
    fprintf(stderr, "    -P FLOAT    Prior probability for heterozygous variant [%1.3f]\n", conf->prior1);
    fprintf(stderr, "    -Q FLOAT    Prior probability for homozygous variant [%1.3f]\n", conf->prior2);
    fprintf(stderr, "    -h          This help\n");
    fprintf(stderr, "\n");

    return 1;
}

int main_pileup(int argc, char *argv[]) {

    int c, i; unsigned j;
    char *reg = 0;
    char *tum = 0;
    char *nor = 0;
    char *outfn = 0;
    char *statsfn = 0;
    pileup_conf_t conf;
    pileup_conf_init(&conf);

    if (argc<2) return usage(&conf);
    while ((c=getopt(argc, argv, ":o:w:g:@:5:3:b:s:E:M:x:C:P:Q:t:n:m:a:l:T:I:SNrcdupv:h"))>=0) {
        switch (c) {
            case 'g': reg = optarg; break;
            case '@': conf.bt.n_threads = atoi(optarg); break;
            case 's': conf.bt.step = atoi(optarg); break;
            case 'N': conf.comm.is_nome = 1; break;
            case 'S': conf.somatic = 1; break;
            case 'T': tum = optarg; break;
            case 'I': nor = optarg; break;

            case 'o': outfn = optarg; break;
            case 'w': statsfn = strdup(optarg); break;
            case 'v': conf.comm.verbose = atoi(optarg); break;

            case 'b': conf.filt.min_base_qual = atoi(optarg); break;
            case 'm': conf.filt.min_mapq = atoi(optarg); break;
            case 'a': conf.filt.min_score = atoi(optarg); break;
            case 't': conf.filt.max_retention = atoi(optarg); break;
            case 'l': conf.filt.min_read_len = atoi(optarg); break;
            case '5': conf.filt.min_dist_end_5p = atoi(optarg); break;
            case '3': conf.filt.min_dist_end_3p = atoi(optarg); break;
            case 'r': conf.ambi_redist = 0; break;
            case 'c': conf.filt.filter_secondary = 0; break;
            case 'd': conf.filt.filter_doublecnt = 0; break;
            case 'u': conf.filt.filter_duplicate = 0; break;
            case 'p': conf.filt.filter_ppair = 0; break;
            case 'n': conf.filt.max_nm = atoi(optarg); break;

            case 'E': conf.error = atof(optarg); break;
            case 'M': conf.mu = atof(optarg); break;
            case 'x': conf.mu_somatic = atof(optarg); break;
            case 'C': conf.contam = atof(optarg); break;
            case 'P': conf.prior1 = atof(optarg); break;
            case 'Q': conf.prior2 = atof(optarg); break;

            case 'h': return usage(&conf);
            case ':': usage(&conf); wzfatal("Option needs an argument: -%c\n", optopt); break;
            case '?': usage(&conf); wzfatal("Unrecognized option: -%c\n", optopt); break;
            default:
                      return usage(&conf);
        }
    }

    if (optind + 2 > argc && !conf.somatic) {
        usage(&conf);
        wzfatal("Reference or bam input is missing\n");
    } else if (optind + 1 > argc && conf.somatic) {
        usage(&conf);
        wzfatal("Reference input is missing\n");
    }

    char *reffn = argv[optind++];
    char **in_fns = 0; int n_fns = 0;
    if (conf.somatic) {
        if (!tum) {
            fprintf(stderr, "[%s:%d] To call somatic events (-S), we need to specify the tumor BAM (-T).\nAbort.\n", __func__, __LINE__);
            fflush(stderr);
            exit(1);
        } else if (!nor) {
            fprintf(stderr, "[%s:%d] To call somatic events (-S), we need to specify the normal BAM (-I).\nAbort.\n", __func__, __LINE__);
            fflush(stderr);
            exit(1);
        }

        n_fns = 2;
        in_fns = realloc(in_fns, (n_fns)*sizeof(char*));
        in_fns[0] = tum;
        in_fns[1] = nor;
    } else {
        if (tum) {
            fprintf(stderr, "[%s:%d] You have specified a tumor BAM (-T) for somatic mode, but -S was not supplied.\nAbort.\n", __func__, __LINE__);
            fflush(stderr);
            exit(1);
        } else if (nor) {
            fprintf(stderr, "[%s:%d] You have specified a normal BAM (-I) for somatic mode, but -S was not supplied.\nAbort.\n", __func__, __LINE__);
            fflush(stderr);
            exit(1);
        }

        for (; optind < argc; ++optind) {
            in_fns = realloc(in_fns, (n_fns+1)*sizeof(char*));
            in_fns[n_fns++] = argv[optind];
        }
    }

    if (conf.comm.verbose > 5) {
        for (i=0; i<n_fns; ++i)
            fprintf(stderr, "[%s:%d] in bams: %s\n", __func__, __LINE__, in_fns[i]);
        fflush(stderr);
    }

    /* read header in the 1st bam, assume all bams are equal in header */
    htsFile *in = hts_open(in_fns[0], "rb");
    if (!in) {
        fprintf(stderr, "[%s:%d] Cannot open %s\nAbort.\n", __func__, __LINE__, in_fns[0]);
        fflush(stderr);
        exit(1);
    }

    // sort sequence name by alphabetic order, chr1, chr10, chr11
    bam_hdr_t *hdr = sam_hdr_read(in);
    target_v *targets = init_target_v(50);
    target_t *t;
    for (i=0; i<hdr->n_targets; ++i) {
        t = next_ref_target_v(targets);
        t->tid = i;
        t->name = hdr->target_name[i];
        t->len = hdr->target_len[i];
    }
    qsort(targets->buffer, targets->size, sizeof(target_t), compare_targets);

    char *vcf_hdr = print_vcf_header(reffn, targets, argv, argc, &conf, in_fns, n_fns);

    // setup writer
    pthread_t writer;
    writer_conf_t writer_conf = {
        .q = wqueue_init(record, 100000),
        .bam_fns = in_fns,
        .n_bams = n_fns,
        .outfn = outfn,
        .statsfn = statsfn,
        .header = vcf_hdr,
        .targets = targets,
        .conf = &conf,
    };
    pthread_create(&writer, NULL, write_func, &writer_conf);

    // send out work
    wqueue_t(window) *wq = wqueue_init(window, 100000);
    pthread_t *processors = calloc(conf.bt.n_threads, sizeof(pthread_t));
    result_t *results = calloc(conf.bt.n_threads, sizeof(result_t));
    for (i=0; i<conf.bt.n_threads; ++i) {
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        results[i].ref_fn = reffn;
        results[i].bam_fns = in_fns;
        results[i].n_bams = n_fns;
        results[i].conf = &conf;
        pthread_create(&processors[i], NULL, process_func, &results[i]);
    }

    window_t w; memset(&w, 0, sizeof(window_t));
    uint32_t wbeg;
    int64_t block_id=0;

    /* process bam */
    if (reg) {                    /* regional */
        int tid;
        uint32_t beg, end;
        biscuit_parse_region(reg, hdr, &tid, (int*) &beg, (int*) &end);
        /* chromosome are assumed to be less than 2**29 */
        beg++; // shift beg from 0-based to 1-based
        if (beg<=0) beg = 1;
        if (end>hdr->target_len[tid]) end = hdr->target_len[tid];
        for (wbeg = beg; wbeg < end; wbeg += conf.bt.step, block_id++) {
            w.tid = tid;
            w.block_id = block_id;
            w.beg = wbeg;
            w.end = wbeg + conf.bt.step;
            if (w.end > end) w.end = end;
            wqueue_put(window, wq, &w);
        }
    } else {                      /* entire bam */
        for (j=0; j<targets->size; ++j) {
            t = ref_target_v(targets, j);
            for (wbeg = 1; wbeg < t->len; wbeg += conf.bt.step, block_id++) {
                w.tid = t->tid;
                w.block_id = block_id;
                w.beg = wbeg;
                w.end = wbeg+conf.bt.step;
                if (w.end > t->len) w.end = t->len;
                wqueue_put(window, wq, &w);
            }
        }
    }
    for (i=0; i<conf.bt.n_threads; ++i) {
        w.tid = -1;
        wqueue_put(window, wq, &w);
    }

    for (i=0; i<conf.bt.n_threads; ++i) {
        pthread_join(processors[i], NULL);
    }

    record_t rec = { .block_id = RECORD_QUEUE_END };
    wqueue_put2(record, writer_conf.q, rec);
    pthread_join(writer, NULL);
    wqueue_destroy(record, writer_conf.q);

    free_target_v(targets);
    free(results);
    free(processors);
    free(vcf_hdr);
    wqueue_destroy(window, wq);
    hts_close(in);
    bam_hdr_destroy(hdr);
    free(in_fns);

    return 0;
}

