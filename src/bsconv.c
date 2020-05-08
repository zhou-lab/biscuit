/* annotate bisulfite conversion of reads
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
#include "encode.h"

typedef struct {
    int max_cph;
    float max_cph_frac;
    int max_cpa;
    int max_cpc;
    int max_cpt;
    int filter_u;
    int show_filtered;
    int print_in_tab;
} bsconv_conf_t;

typedef struct bsconv_data_t {
    refcache_t *rs;
    bsconv_conf_t *conf;
    int n;
    int n_filtered;
} bsconv_data_t;

static int bsconv_func(bam1_t *b, samFile *out, bam_hdr_t *hdr, void *data) {

    (void) (out);
    bsconv_data_t *d = (bsconv_data_t*) data;
    bsconv_conf_t *conf = d->conf;
    const bam1_core_t *c = &b->core;

    int retn[5] = {0};
    int conv[5] = {0};

    int tofilter = 0;
    if (c->flag & BAM_FUNMAP ||
            c->flag & BAM_FQCFAIL)
        /* c->flag & BAM_FSECONDARY || */
        /* c->flag & BAM_FDUP || */
        /* c->flag & BAM_FSUPPLEMENTARY */
    {
        tofilter = 1;
        goto OUTPUT;
    }

    // TODO: this requires "-" input be input with "samtools view -h", drop this
    refcache_fetch(d->rs, hdr->target_name[c->tid], max(1,c->pos-10), bam_endpos(b)+10);
    uint32_t rpos=c->pos+1, qpos=0;
    int i; unsigned j;
    char rb, qb;
    uint8_t bsstrand = get_bsstrand(d->rs, b, 0, conf->filter_u);
    if (bsstrand == 2) {
        tofilter = 1;
        goto OUTPUT;
    }
    char fivenuc[5];

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

                    // count retention and conversion
                    qb = toupper(bscall(b, qpos+j));
                    if (bsstrand && rb == 'G') {
                        if (qb == 'G') retn[nt256char_to_nt256int8_table[(unsigned char)fivenuc[3]]]++;
                        else if (qb == 'A') conv[nt256char_to_nt256int8_table[(unsigned char)fivenuc[3]]]++;
                    } else if (!bsstrand && rb == 'C') {
                        if (qb == 'C') retn[nt256char_to_nt256int8_table[(unsigned char)fivenuc[3]]]++;
                        else if (qb == 'T') conv[nt256char_to_nt256int8_table[(unsigned char)fivenuc[3]]]++;
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
            case BAM_CHARD_CLIP:
                qpos += oplen;
                break;
            default:
                fprintf(stderr, "Unknown cigar, %u\n", op);
                abort();
        }
    }

    if (conf->max_cpa >= 0 && retn[nt256char_to_nt256int8_table['A']] > conf->max_cpa) tofilter = 1;
    if (conf->max_cpc >= 0 && retn[nt256char_to_nt256int8_table['C']] > conf->max_cpc) tofilter = 1;
    if (conf->max_cpt >= 0 && retn[nt256char_to_nt256int8_table['T']] > conf->max_cpt) tofilter = 1;
    if (conf->max_cph >= 0 &&
            retn[nt256char_to_nt256int8_table['A']] +
            retn[nt256char_to_nt256int8_table['C']] +
            retn[nt256char_to_nt256int8_table['T']] > conf->max_cph) tofilter = 1;

    if (conf->max_cph_frac < 1.0) {
        int cph_retn = retn[nt256char_to_nt256int8_table['A']] +
            retn[nt256char_to_nt256int8_table['C']] +
            retn[nt256char_to_nt256int8_table['T']];
        int cph_conv = conv[nt256char_to_nt256int8_table['A']] +
            conv[nt256char_to_nt256int8_table['C']] +
            conv[nt256char_to_nt256int8_table['T']];
        if (cph_retn + cph_conv > 0 &&
                (float) cph_retn / (cph_retn + cph_conv) > conf->max_cph_frac) tofilter = 1;
    }

OUTPUT:
    d->n++;
    if (tofilter) d->n_filtered++;
    if (conf->show_filtered) tofilter = !tofilter;
    if (tofilter) return 0;

    // tab output
    if (conf->print_in_tab) {
        for (i=0; i<4; ++i) {
            if (i) putchar('\t');
            printf("%d\t%d", retn[i], conv[i]);
        }
        printf("\t%s", bam_get_qname(b));
        putchar('\n');
        return 0;
    }

    // bam output
    kstring_t s;
    if (!conf->print_in_tab && out) {
        s.m = s.l = 0; s.s = 0;
        for (i=0; i<4; ++i) {
            if (i) kputc(',', &s);
            ksprintf(&s, "C%c_R%dC%d", nt256int8_to_nt256char_table[i], retn[i], conv[i]);
        }
        bam_aux_append(b, "ZN", 'Z', s.l+1, (uint8_t*) s.s);
        free(s.s);

        if (sam_write1(out, hdr, b) < 0) {
            wzfatal("Cannot write bam.\n");
        }
    } else {
        s.m = s.l = 0; s.s = 0;
        sam_format1(hdr, b, &s);
        if (puts(s.s) < 0)
            if (errno == EPIPE) exit(1);
        free(s.s);
    }
    return 0;
}

static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit bsconv [options] <ref.fa> <in.bam> [out.bam]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -g STR      Region (optional, will process the whole bam if not specified)\n");
    fprintf(stderr, "    -u          Filter unclear bs-strand (YD:u) reads\n");
    fprintf(stderr, "    -m INT      Filter: maximum CpH retention [Inf]\n");
    fprintf(stderr, "    -f FLOAT    Filter: maximum CpH retention fraction [1.0]\n");
    fprintf(stderr, "    -a INT      Filter: maximum CpA retention [Inf]\n");
    fprintf(stderr, "    -c INT      Filter: maximum CpC retention [Inf]\n");
    fprintf(stderr, "    -t INT      Filter: maximum CpT retention [Inf]\n");
    fprintf(stderr, "    -p          Print in tab-separated format, print order:\n");
    fprintf(stderr, "                    CpA_R, CpA_C, CpC_R, CpC_C, CpG_R, CpG_C, CpT_R, CpT_C\n");
    fprintf(stderr, "    -v          Show filtered reads instead of remaining reads\n");
    fprintf(stderr, "    -h          This help\n");
    fprintf(stderr, "\n");
}

int main_bsconv(int argc, char *argv[]) {
    int c;
    char *reg = 0; // target region
    bsconv_conf_t conf = {0};
    conf.max_cph = conf.max_cpa = conf.max_cpc = conf.max_cpt = -1;
    conf.max_cph_frac = 1.0;
    conf.print_in_tab = 0;

    if (argc < 2) { usage(); return 1; }
    while ((c = getopt(argc, argv, ":g:m:ac:f:pt:uvh")) >= 0) {
        switch (c) {
            case 'g': reg = optarg; break;
            case 'm': conf.max_cph = atoi(optarg); break;
            case 'f': conf.max_cph_frac = atof(optarg); break;
            case 'a': conf.max_cpa = atoi(optarg); break;
            case 'c': conf.max_cpc = atoi(optarg); break;
            case 't': conf.max_cpt = atoi(optarg); break;
            case 'u': conf.filter_u = 1; break;
            case 'p': conf.print_in_tab = 1; break;
            case 'v': conf.show_filtered = 1; break;
            case 'h': usage(); return 1;
            case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt);
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt);
            default: usage(); return 1;
        }
    }

    char *reffn = optind < argc ? argv[optind++] : NULL;
    char *infn  = optind < argc ? argv[optind++] : NULL;
    char *outfn = optind < argc ? argv[optind++] : NULL;
    if (!reffn || !infn) {
        usage();
        wzfatal("Please provide reference and input bam.\n");
    }

    bsconv_data_t d = {0};
    d.rs = init_refcache(reffn, 100, 100000);
    d.conf = &conf;
    bam_filter(infn, outfn, reg, &d, bsconv_func);

    fprintf(stderr, "\n");
    fprintf(stderr, "[%s:%d] Processed %d reads, %d (%f%%) remains.\n",
            __func__, __LINE__, d.n, d.n - d.n_filtered, (double) (d.n - d.n_filtered) / d.n * 100);
    fflush(stderr);

    free_refcache(d.rs);
    return 0;
}
