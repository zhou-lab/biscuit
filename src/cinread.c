/* Print C in reads as long form
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Wanding.Zhou@vai.org
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

#include <unistd.h>
#include <errno.h>
#include "wstr.h"
#include "wzmisc.h"
#include "refcache.h"
#include "sam.h"
#include "bamfilter.h"
#include "pileup.h"

const char *tp_names[] = {
    "QNAME", "QPAIR", "STRAND", "BSSTRAND", "MAPQ",
    "QBEG", "QEND",
    "CHRM",
    "CRPOS",
    "CGRPOS",
    "CQPOS",
    "CRBASE",
    "CCTXT",
    "CQBASE",
    "CRETENTION"};

typedef enum {
    TP_QNAME, TP_QPAIR, TP_STRAND, TP_BSSTRAND, TP_MAPQ,
    TP_QBEG, TP_QEND, // query start, end
    TP_CHRM,          // chromosome
    TP_CRPOS,         // cytosine position on reference
    TP_CGRPOS,        // CpG position on reference (-1 if not applicable)
    TP_CQPOS,         // cytosine position on query
    TP_CRBASE,        // cytosine reference
    TP_CCTXT,         // cytosine context, strand flipped
    TP_CQBASE,        // base called on read
    TP_CRETENTION} __tp_name_t; // retention (R) or conversion (C)

const char *tgt_names[] = {"c", "cg", "ch", "hcg", "gch", "hch"};

typedef enum {SL_C, SL_CG, SL_CH, SL_HCG, SL_GCH, SL_HCH} __tgt_name_t; // SL_ select target

typedef struct {
    int n_tp_names;
    __tp_name_t *tp_names;
    __tgt_name_t tgt;
    FILE *out;
    int skip_secondary;
} cinread_conf_t;

typedef struct cinread_data_t {
    refcache_t *rs;
    FILE *output;
    cinread_conf_t *conf;
} cinread_data_t;

static int cinread_func(bam1_t *b, samFile *out, bam_hdr_t *hdr, void *data) {

    (void) (out);
    cinread_data_t *d = (cinread_data_t*) data;
    cinread_conf_t *conf = d->conf;
    const bam1_core_t *c = &b->core;

    if (c->flag & BAM_FUNMAP) return 0; // skip unmapped
    if (conf->skip_secondary && c->flag & BAM_FSECONDARY) return 0; // skip secondary

    // TODO: this requires "-" input be input with "samtools view -h", drop this
    refcache_fetch(d->rs, hdr->target_name[c->tid], max(1,c->pos-10), bam_endpos(b)+10);
    uint32_t rpos=c->pos+1, qpos=0;
    int i, k; unsigned j;
    char rb, qb;
    char retention = 'N';
    uint8_t bsstrand = get_bsstrand(d->rs, b, 0, 0);
    char fivenuc[5];
    int is_tgt = 0;
    int l_qseq = c->l_qseq;
    for (i=0; i<c->n_cigar; ++i) {
        uint32_t op = bam_cigar_op(bam_get_cigar(b)[i]);
        uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
        switch(op) {
            case BAM_CMATCH:
                for(j=0; j<oplen; ++j) {
                    rb = refcache_getbase_upcase(d->rs, rpos+j);

                    // avoid looking at the wrong strand
                    if (rb != 'C' && rb != 'G') continue;
                    if (bsstrand && rb == 'C') continue;
                    if (!bsstrand && rb == 'G') continue;

                    // filter by context
                    fivenuc_context(d->rs, rpos+j, rb, fivenuc);
                    is_tgt = 0;
                    switch (conf->tgt) {
                        case SL_C: is_tgt=1; break;
                        case SL_CG: if (fivenuc[3] == 'G') is_tgt=1; break;
                        case SL_CH: if (fivenuc[3] != 'G') is_tgt=1; break;
                        case SL_HCG: if (fivenuc[3] == 'G' && fivenuc[1] != 'G') is_tgt=1; break;
                        case SL_GCH: if (fivenuc[3] != 'G' && fivenuc[1] == 'G') is_tgt=1; break;
                        case SL_HCH: if (fivenuc[3] != 'G' && fivenuc[1] != 'G') is_tgt=1; break;
                        default: wzfatal("Unknown target name: %u\n", conf->tgt);
                    }
                    if (!is_tgt) continue;

                    // set retention
                    qb = toupper(bscall(b, qpos+j));
                    if (bsstrand && rb == 'G') {
                        if (qb == 'G') retention = 'R';
                        else if (qb == 'A') retention = 'C';
                        else retention = 'N';
                    } else if (!bsstrand && rb == 'C') {
                        if (qb == 'C') retention = 'R';
                        else if (qb == 'T') retention = 'C';
                        else retention = 'N';
                    } else retention = 'N';

                    for (k=0; k<conf->n_tp_names; ++k) {
                        if (k) fputc('\t', conf->out);
                        switch(conf->tp_names[k]) {
                            case TP_QNAME: fputs(bam_get_qname(b), conf->out); break;
                            case TP_QPAIR: fputc((c->flag&BAM_FREAD2)?'2':'1', conf->out); break;
                            case TP_QBEG: fprintf(conf->out, "%d", c->pos+1); break;
                            case TP_QEND: fprintf(conf->out, "%d", bam_endpos(b)); break;
                            case TP_STRAND: fputc((c->flag&BAM_FREVERSE)?'-':'+', conf->out); break;
                            case TP_BSSTRAND: fputc(bsstrand?'-':'+', conf->out); break;
                            case TP_MAPQ: fprintf(conf->out, "%d", c->qual); break;
                            case TP_CHRM: fputs(hdr->target_name[c->tid], conf->out); break;
                            case TP_CRPOS: fprintf(conf->out, "%u", rpos+j); break;
                            case TP_CGRPOS: {
                                if (fivenuc[3] == 'G') {
                                    if (rb == 'C') fprintf(conf->out, "%u", rpos+j);
                                    else if (rb == 'G') fprintf(conf->out, "%u", rpos+j-1);
                                } else fprintf(conf->out, "-1");
                                break;
                            }
                            case TP_CQPOS: {
                                // note when there is hard clipping, l_qseq might be < qpos+j
                                // see following for compensation
                                fprintf(conf->out, "%u", (c->flag&BAM_FREVERSE)?(l_qseq-qpos-j):(qpos+j));
                                break;
                            }
                            case TP_CRBASE: fputc(rb, conf->out); break;
                            case TP_CCTXT: {
                                fprintf(conf->out, "%.5s", fivenuc);
                                break;
                            }
                            case TP_CQBASE: fputc(qb, conf->out); break;
                            case TP_CRETENTION: fputc(retention, conf->out); break;
                            default: wzfatal("Unknown print name: %u\n", conf->tp_names[k]);
                        }
                    }

                    if (fputc('\n', conf->out) < 0)
                        if (errno == EPIPE) exit(1);
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
                l_qseq += oplen; // c->l_qseq excludes hard clipping, add back here.
                break;
            default:
                fprintf(stderr, "Unknown cigar, %u\n", op);
                abort();
        }
    }

    return 0;
}


static void usage() {

    unsigned i;
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: cinread [options] ref.fa in.bam\n");
    fprintf(stderr, "Input options:\n");
    fprintf(stderr, "     -g        region.\n");
    fprintf(stderr, "     -t        target (");
    for (i=0; i<sizeof(tgt_names)/sizeof(tgt_names[0]); ++i) {
        if (i) fputs(", ", stderr);
        fputs(tgt_names[i], stderr);
    } fputs(") [cg]\n", stderr);
    fprintf(stderr, "     -p        content to print, \",\"-delimited:\n");
    for(i=0; i<sizeof(tp_names)/sizeof(tp_names[0]); ++i) {
        if (i%5 == 0) fputs("\n               ", stderr);
        else fputs(", ", stderr);
        fputs(tp_names[i], stderr);
    } fputs("\n\n", stderr);
    fputs("               [QNAME,QPAIR,BSSTRAND,CRBASE,CQBASE]\n\n", stderr);
    fprintf(stderr, "     -s        consider secondary mapping.\n");
    fprintf(stderr, "     -o        output.\n");
    fprintf(stderr, "     -h        this help.\n");
    fprintf(stderr, "\n");
}

int main_cinread(int argc, char *argv[]) {
    int c;
    char *reg = 0; // target region
    cinread_conf_t conf = {0};
    conf.skip_secondary = 1;
    char *outfn = NULL;

    char *tgt_str = 0; char *tp_str = 0;
    if (argc < 2) { usage(); return 1; }
    while ((c = getopt(argc, argv, "g:o:t:p:sh")) >= 0) {
        switch (c) {
            case 'g': reg = optarg; break;
            case 'o': outfn = optarg; break;
            case 't': tgt_str = optarg; break;
            case 'p': tp_str = optarg; break;
            case 's': conf.skip_secondary = 0; break;
            case 'h': usage(); return 1;
            default:
                fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
                exit(1);
                break;
        }
    }

    // parse target base
    unsigned i;
    if (tgt_str) {
        for (i=0; i<sizeof(tgt_names)/sizeof(tgt_names[0]); ++i) {
            if (strcmp(tgt_names[i], tgt_str)==0) {
                conf.tgt = i; break;
            }
        }
        if (i == sizeof(tgt_names)/sizeof(tgt_names[0])) {
            usage();
            wzfatal("Target name %s unrecognized.\n", tgt_str);
        }
    } else {
        conf.tgt = SL_CG;
    }

    // parse column to print
    if (tp_str) {
        char *pch = strtok(tp_str, ",");
        while (pch != NULL) {
            for (i=0; i<sizeof(tp_names)/sizeof(tp_names[0]); ++i) {
                if (strcmp(pch, tp_names[i]) == 0) {
                    conf.tp_names = realloc(conf.tp_names, (conf.n_tp_names+1)*sizeof(__tp_name_t));
                    conf.tp_names[conf.n_tp_names++] = i;
                    break;
                }
            }
            if (i == sizeof(tp_names)/sizeof(tp_names[0])) {
                usage();
                wzfatal("Print name %s unrecognized.\n", pch);
            }
            pch = strtok(NULL, ",");
        }
    } else {
        conf.n_tp_names = 5;
        conf.tp_names = malloc((conf.n_tp_names)*sizeof(__tp_name_t));
        conf.tp_names[0] = TP_QNAME;
        conf.tp_names[1] = TP_QPAIR;
        conf.tp_names[2] = TP_BSSTRAND;
        conf.tp_names[3] = TP_CRBASE;
        conf.tp_names[4] = TP_CQBASE;
    }


    /* fprintf(stderr, "conf target: %d - %s\n", conf.tgt, tgt_names[conf.tgt]); */
    /* for (i=0; i<conf.n_tp_names; ++i) */
    /*   fprintf(stderr, "conf to print: %d - %s\n", i, tp_names[conf.tp_names[i]]); */

    char *reffn = optind < argc ? argv[optind++] : NULL;
    char *infn  = optind < argc ? argv[optind++] : NULL;
    if (!reffn || !infn) {
        usage();
        wzfatal("Please provide reference and input bam.\n");
    }

    if (outfn) conf.out = fopen(outfn, "w");
    else conf.out = stdout;

    cinread_data_t d = {0};
    d.rs = init_refcache(reffn, 100, 100000);
    d.conf = &conf;
    bam_filter(infn, 0, reg, &d, cinread_func);

    free_refcache(d.rs);
    free(conf.tp_names);
    if (outfn) fclose(conf.out);
    return 0;
}
