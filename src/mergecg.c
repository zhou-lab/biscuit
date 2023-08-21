/* Merge C and G beta values in CpG dinucleotide context
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021-2023 Jacob.Morrison@vai.org
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
**/

#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <string.h>
#include <zlib.h>
#include "wzbed.h"
#include "refcache.h"

typedef struct bed_data_meth_t {
    char ref;
    int nsamples;
    double *c_betas;
    double *g_betas;
    int *c_depts;
    int *g_depts;
} bed_data_meth_t;

void init_data_meth(bed1_t *b, void *aux_data) {
    (void) aux_data;
    b->data = calloc(1, sizeof(bed_data_meth_t));
}

void free_data_meth(void *bed_data) {
    free(((bed_data_meth_t*) bed_data)->c_betas);
    free(((bed_data_meth_t*) bed_data)->c_depts);
    free(((bed_data_meth_t*) bed_data)->g_betas);
    free(((bed_data_meth_t*) bed_data)->g_depts);
    free(bed_data);
}

void parse_data_meth(bed1_t *b, char **fields, int nfields) {

    int start;
    if (strcmp(fields[3],"C")==0 || strcmp(fields[3],"G")==0) start = 7;
    else start = 3;

    bed_data_meth_t *bd = (bed_data_meth_t*) b->data;
    if (bd->nsamples <= 0) { // first line
        bd->nsamples = (nfields-start)/2;
        if (bd->nsamples <= 0) wzfatal("No sample data identified.\n");
        bd->c_betas = calloc(bd->nsamples, sizeof(double));
        bd->c_depts = calloc(bd->nsamples, sizeof(int));
        bd->g_betas = calloc(bd->nsamples, sizeof(double));
        bd->g_depts = calloc(bd->nsamples, sizeof(int));
    } else if (bd->nsamples*2 + start != nfields) wzfatal("Malformed bed file.\n");

    // assume it's a C
    int i;
    for (i=0; i<bd->nsamples; ++i) {
        bd->c_betas[i] = atof(fields[start+2*i]);
        bd->c_depts[i] = atoi(fields[start+1+2*i]);
        memset(bd->g_betas, 0, bd->nsamples*sizeof(double));
        memset(bd->g_depts, 0, bd->nsamples*sizeof(int));
    }
}

void format_output(bed1_t *p, char *chrm, char base_before, char base_after, int min_depth) {
    // skip reporting if all sample has zero coverage
    int i; int max_depth = 0;
    bed_data_meth_t *pd = (bed_data_meth_t*) p->data;
    for (i=0; i<pd->nsamples; ++i) {
        if (pd->c_depts[i] + pd->g_depts[i] > max_depth) {
            max_depth = pd->c_depts[i] + pd->g_depts[i];
        }
    }
    if (max_depth == 0 || max_depth < min_depth) return;

    // if p is in CpG contxt, adjust the coordinates
    if (pd->ref == 'C' && base_after == 'G') p->end++;
    else if (pd->ref == 'G' && base_before == 'C') p->beg--;

    printf("%s\t%"PRId64"\t%"PRId64, chrm, p->beg, p->end);
    for (i=0; i<pd->nsamples; ++i) {
        int cov = pd->c_depts[i] + pd->g_depts[i];
        if (cov == 0) {
            fputs("\t.\t0", stdout);
        } else {
            // Calculate merged beta value from count values to reduce error in beta value
            float c_ret = rintf(pd->c_betas[i]*pd->c_depts[i]);
            float g_ret = rintf(pd->g_betas[i]*pd->g_depts[i]);
            printf("\t%1.3f\t%d", (c_ret + g_ret) / (double) cov, cov);
        }
        if (pd->c_depts[i] == 0) {
            fputs("\tC:.:0", stdout);
        } else {
            printf("\tC:%1.3f:%d", pd->c_betas[i], pd->c_depts[i]);
        }
        if (pd->g_depts[i] == 0) {
            fputs(",G:.:0", stdout);
        } else {
            printf(",G:%1.3f:%d", pd->g_betas[i], pd->g_depts[i]);
        }
    }
    putchar('\n');
}

static int usage(int min_depth) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit mergecg [options] <ref.fa> <in.bed>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -N        NOMe-seq mode, only merge C,G both in HCGD context\n");
    fprintf(stderr, "    -k INT    Minimum depth after merging - applies to the maximum depth\n");
    fprintf(stderr, "                  across samples [%d]\n", min_depth);
    fprintf(stderr, "    -h        This help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note, in.bed is a position sorted bed file with beta values and coverages found\n");
    fprintf(stderr, "    in columns 4 and 5, respectively. Additional beta value-coverage column\n");
    fprintf(stderr, "    pairs are added for each additional sample. This is the format that would be\n");
    fprintf(stderr, "    found in the output of biscuit vcf2bed without the '-e' flag included.\n");
    fprintf(stderr, "\n");

    return 1;
}

int main_mergecg(int argc, char *argv[]) {

    int c;
    int nome_mode = 0;
    int min_depth = 0;
    if (argc<2) return usage(min_depth); 
    while ((c = getopt(argc, argv, ":k:hN"))>=0) {
        switch (c) {
            case 'N': nome_mode=1; break;
            case 'k': min_depth = atoi(optarg); break;
            case 'h': return usage(min_depth); break;
            case ':': usage(min_depth); wzfatal("Option needs an argument: -%c\n", optopt); break;
            case '?': usage(min_depth); wzfatal("Unrecognized option: -%c\n", optopt); break;
            default: usage(min_depth); wzfatal("Unrecognized option: %c.\n", c); break;
        }
    }

    if (optind + 2 > argc) { 
        usage(min_depth); 
        wzfatal("Please supply reference file and sorted bed file.\n"); 
    }

    refcache_t *rc = init_refcache(argv[optind++], 10000, 10000);
    bed_file_t *bed = init_bed_file(argv[optind++]);
    bed1_t *b = init_bed1(init_data_meth, NULL);
    bed1_t *p = init_bed1(init_data_meth, NULL);
    b->tid = -1; p->tid = -1;
    bed_data_meth_t *bd, *pd;
    char p_base_before='N', p_base_after='N';
    char b_base_before='N', b_base_after='N';
    while (bed_read1(bed, b, parse_data_meth)) {
        refcache_fetch_chrm(rc, tid2name(bed->targets, b->tid));
        bd = (bed_data_meth_t*) b->data;
        bd->ref = refcache_getbase_upcase(rc, b->end);
        if (b->end-1 < 0) b_base_before = 'N';
        else b_base_before = refcache_getbase_upcase(rc, b->end-1);
        if (b->end == rc->end) b_base_after = 'N';
        else b_base_after = refcache_getbase_upcase(rc, b->end+1);
        if (bd->ref == 'G') { // correct based on the actual reference base
            memcpy(bd->g_betas, bd->c_betas, bd->nsamples*sizeof(double));
            memcpy(bd->g_depts, bd->c_depts, bd->nsamples*sizeof(int));
            memset(bd->c_betas, 0, bd->nsamples*sizeof(double));
            memset(bd->c_depts, 0, bd->nsamples*sizeof(int));
        }

        // merge C, G
        pd = (bed_data_meth_t*) p->data;
        if (b->tid == p->tid && 
                b->beg == p->beg+1 && b->end == p->end+1 && 
                bd->ref == 'G' && pd->ref == 'C' &&
                (!nome_mode ||
                 (p_base_before != 'G' && b_base_after != 'C'))) {
            if (pd->nsamples != bd->nsamples) wzfatal("Missing sample at %s:%u-%u.\n", rc->chrm, b->beg, b->end);
            memcpy(pd->g_betas, bd->g_betas, bd->nsamples*sizeof(double));
            memcpy(pd->g_depts, bd->g_depts, bd->nsamples*sizeof(int));
            b->tid = -1; // merged, no more reporting
        }

        // report p
        if (p->tid >= 0) format_output(p, tid2name(bed->targets, p->tid), p_base_before, p_base_after, min_depth);
        bed1_t *tmp = p; p = b; b = tmp;
        p_base_before = b_base_before; p_base_after = b_base_after;
    }

    // report last p
    if (p->tid >= 0) format_output(p, tid2name(bed->targets, p->tid), p_base_before, p_base_after, min_depth);
    free_bed1(b, free_data_meth);
    free_bed1(p, free_data_meth);
    free_refcache(rc);
    free_bed_file(bed);

    return 0;
}
