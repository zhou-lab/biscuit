/* Rectangularize an epiread output 
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
 * 
 **/

#include "wztsv.h"
#include "wvec.h"
#include "kstring.h"
#include "refcache.h"

int refcache_next_cg(refcache_t *rc, int pos) {
    for (;; pos++) {
        if (toupper(refcache_getbase_auto(rc, NULL, pos)) == 'C' &&
                toupper(refcache_getbase_auto(rc, NULL, pos+1)) == 'G')
            return pos;
    }
}

typedef struct read_t {
    kstring_t seq;
    kstring_t other; // other fields of the line
} read_t;

DEFINE_VECTOR(read_v, read_t)

static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit rectangle [options] <ref.fa> <in.epiread>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR    Output file [stdout]\n");
    fprintf(stderr, "    -h        This help\n");
    fprintf(stderr, "Note, this is not currently compatible with epiread format when run with the -A flag\n");
    fprintf(stderr, "    i.e., biscuit epiread -A [-B snps.bed] <ref.fa> <in.bam>\n");
    fprintf(stderr, "\n");
}

int main_rectangle(int argc, char *argv[]) {

    char *out_fn = 0;

    if (argc<2) { usage(); return 1; }
    int c;
    while ((c = getopt(argc, argv, ":o:h"))>=0) {
        switch (c) {
            case 'o': out_fn = optarg; break;
            case 'h': usage(); return 1;
            case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt); break;
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt); break;
            default: usage(); return 1;
        }
    }

    if (optind + 2 > argc) {
        usage();
        wzfatal("Reference file or epiread file is missing\n");
    }
    char *ref_fn = argv[optind++];
    char *epiread_fn = argv[optind++];
    refcache_t *rc = init_refcache(ref_fn, 1000, 1000);

    tsv_t *tsv = tsv_open(epiread_fn);
    int region_beg = 0, region_width = -1;
    read_v *reads = init_read_v(1000);
    char *chrm = NULL;
    while (tsv_read(tsv)) {
        if (tsv_is_blankline(tsv)) continue;
        char *f = tsv_field(tsv, 4);
        if (f[0] == '.') {
            read_t *r = next_ref_read_v(reads);
            r->seq = (const kstring_t) {0};
            r->other = (const kstring_t) {0};
            kputs(tsv->line, &r->other);
            continue;
        }

        int read_beg = atoi(f);
        if (!region_beg) region_beg = read_beg;

        if (chrm == NULL) {
            chrm = calloc(strlen(tsv_field(tsv, 0))+1, sizeof(char));
            strcpy(chrm, tsv_field(tsv, 0));
            refcache_set_chromosome(rc, chrm);
        } else if (strcmp(chrm, tsv_field(tsv, 0)) != 0) {
            fprintf(
                    stderr,
                    "[%s:%d] Error, rectangle cannot cross chromosomes.\n",
                    __func__, __LINE__);
            exit(1);
        }

        // compute padding
        int p, pad = 0;
        for (p = region_beg; p < read_beg; ++pad)
            p = refcache_next_cg(rc, p) + 1;

        read_t *r = next_ref_read_v(reads);
        r->seq = (const kstring_t) {0};
        while(pad--) kputc('N', &r->seq);
        kputs(tsv->fields[5], &r->seq);
        r->other = (const kstring_t) {0};
        kputs(tsv->line, &r->other);

        // update region width
        if (region_width < 0 || (unsigned) region_width < r->seq.l)
            region_width = r->seq.l;
    }
    free(chrm);

    FILE *out = wzopen_out(out_fn);
    unsigned i;
    for (i = 0; i < reads->size; ++i) {
        read_t *r = ref_read_v(reads, i);
        while (r->seq.l < (unsigned) region_width) kputc('N', &r->seq);
        fputs(r->other.s, out);
        fputc('\t', out);
        fputs(r->seq.s, out);
        fputc('\n', out);
        free(r->seq.s);
        free(r->other.s);
    }

    tsv_close(tsv);
    free_read_v(reads);
    free_refcache(rc);

    return 0;
}
