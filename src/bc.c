/* extract barcodes from FASTQ file
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2023 Jacob.Morrison@vai.org
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

#include "bc.h"

static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit bc [options] <FASTQ 1> <FASTQ 2>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -s    Run for single-end data\n");
    fprintf(stderr, "    -h    This help\n");
    fprintf(stderr, "\n");
}

int main_bc(int argc, char *argv[]) {
    int c;
    bc_conf_t conf = {0};
    conf.single_end = 0;

    if (argc < 2) { usage(); return 1; }
    while ((c = getopt(argc, argv, ":hs")) >= 0) {
        switch(c) {
            case 's':
                conf.single_end = 1;
                break;
            case 'h':
                usage();
                return 0;
            case ':':
                usage();
                fprintf(stderr, "Option needs an argument: -%c\n", optopt);
                return 1;
            case '?':
                usage();
                fprintf(stderr, "Unrecognized option: -%c\n", optopt);
                return 1;
            default:
                usage();
                return 1;
        }
    }

    // Init files and handle read errors
    char *fq1, *fq2;
    kseq_t *ks1, *ks2;
    gzFile fh1 = 0, fh2 = 0;
    fq1 = optind < argc ? argv[optind++] : NULL;
    if (!fq1) {
        usage();
        fprintf(stderr, "No read 1 FASTA file provided\n");
        return 1;
    }

    fh1 = gzopen(fq1, "r");
    if (!fh1) {
        fprintf(stderr, "Could not open input file: %s\n", fq1);
        return 1;
    }

    ks1 = kseq_init(fh1);

    if (!conf.single_end) {
        fq2 = optind < argc ? argv[optind++] : NULL;
        if (!fq2) {
            usage();
            fprintf(stderr, "No read 2 FASTA file provided\n");
            return 1;
        }

        fh2 = gzopen(fq2, "r");
        if (!fh2) {
            fprintf(stderr, "Could not open input file: %s\n", fq2);
            return 1;
        }

        ks2 = kseq_init(fh2);
    } else {
        ks2 = NULL;
    }

    if (fh1) {
        gzclose(fh1);
    }
    if (fh2) {
        gzclose(fh2);
    }

    while (kseq_read(ks1) >= 0) {
        if (ks2 && kseq_read(ks2) < 0) { // read 2 has fewer reads
            fprintf(stderr, "read 2 has fewer sequences\n");
            break;
        }

        fprintf(stdout, "%s %s\n", ks1->name.s, ks2->name.s);
    }
    if (ks2 && kseq_read(ks2) >= 0) {
        fprintf(stderr, "read 1 has fewer sequences\n");
    }

    return 0;
}
