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

gzFile setup_output(const char *ofile, uint8_t read_num) {
    char *fname = (char *)calloc(strlen(ofile)+N_EXTRA, sizeof(char));
    memcpy(fname, ofile, strlen(ofile));
    if (read_num == 0) {
        strcat(fname, ".fq.gz");
    } else if (read_num == 1) {
        strcat(fname, "_R1.fq.gz");
    } else {
        strcat(fname, "_R2.fq.gz");
    }

    gzFile out = gzopen(fname, "w6");
    if (!out) {
        fprintf(stderr, "ERROR: could not open output file: %s\n", ofile);
        free(fname);
        return NULL;
    }

    free(fname);

    return out;
}

int prepare_read(kseq_t *k, kstring_t *s, bc_conf_t *conf, uint8_t remove_bc) {
    // Pre-allocate space, or expand space ahead of time, to reduce the number of allocations
    size_t str_len = k->name.l + k->comment.l + k->seq.l + k->qual.l + (size_t)N_EXTRA;

    if (str_len > s->m) {
        int err = ks_resize(s, str_len);
        if (err < 0) {
            fprintf(stderr, "ERROR: Unable to allocate sufficient space\n");
            return -1;
        }
    }

    // Read name
    ksprintf(s, "@%s", k->name.s);

    // Read comment
    ksprintf(s, " %s", k->comment.s);

    // Sequence
    ksprintf(s, "\n%s\n+\n", k->seq.s);

    // Quality
    ksprintf(s, "%s\n", k->qual.s);

    return 0;
}

void extract_barcodes(bc_conf_t *conf, kseq_t *ks1, kseq_t *ks2, gzFile oh1, gzFile oh2) {
    // Allocate strings
    kstring_t *s1 = (kstring_t *)calloc(1, sizeof(kstring_t));
    kstring_t *s2 = (kstring_t *)calloc(1, sizeof(kstring_t));

    // Process reads
    while (kseq_read(ks1) >= 0) {
        if (ks2 && kseq_read(ks2) < 0) { // read 2 has fewer reads
            fprintf(stderr, "WARNING: read 2 has fewer sequences\n");
            break;
        }

        if (conf->bc_start + conf->bc_length > ks1->seq.l) {
            fprintf(stderr, "WARNING: read is too short to extract barcode, dropping\n");
            continue;
        }

        uint8_t r1_has_bc = conf->mate == 1 ? 1 : 0;
        int err = prepare_read(ks1, s1, conf, r1_has_bc);
        if (err < 0) {
            break;
        }
        if (oh1) {
            if (gzwrite(oh1, s1->s, s1->l) == 0) {
                fprintf(stderr, "ERROR: unsuccessful write to FASTQ\n");
                break;
            }
        } else {
            fprintf(stdout, "%s", s1->s);
        }

        // Reset kstring for next read
        s1->s[0] = '\0';
        s1->l = 0;

        if (ks2) {
            int err = prepare_read(ks2, s2, conf, !r1_has_bc);
            if (err < 0) {
                break;
            }
            if (oh2) {
                if (gzwrite(oh2, s2->s, s2->l) == 0) {
                    fprintf(stderr, "ERROR: unsuccessful write to FASTQ\n");
                    break;
                }
            } else {
                fprintf(stdout, "%s", s2->s);
            }

            // Reset kstring for next read
            s2->s[0] = '\0';
            s2->l = 0;
        }
    }
    if (ks2 && kseq_read(ks2) >= 0) {
        fprintf(stderr, "WARNING: read 1 has fewer sequences\n");
    }

    // Clean up
    free(s2->s);
    free(s2);
    free(s1->s);
    free(s1);
}

static void usage() {
    bc_conf_t conf = {0};
    bc_conf_init(&conf);

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit bc [options] <FASTQ 1> [FASTQ 2]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Processing Options:\n");
    fprintf(stderr, "    -m, --mate INT             which mate the barcode is in (1 or 2) [%u]\n", conf.mate);
    fprintf(stderr, "    -s, --bc-start INT         start position of barcode in read (1-based) [%u]\n", conf.bc_start);
    fprintf(stderr, "    -l, --bc-length INT        length of barcode [%u]\n", conf.bc_length);
    fprintf(stderr, "Output Options:\n");
    fprintf(stderr, "    -o, --output-prefix STR    prefix for output files (NULL writes to stdout) [NULL]\n");
    fprintf(stderr, "    -z, --gz                   write output files as gzipped FASTQs [off]\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "    -h, --help             This help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: When writing to stdout, reads 1 and 2 will alternate (i.e., are interleaved)\n");
    fprintf(stderr, "\n");
}

int main_bc(int argc, char *argv[]) {
    int c;
    bc_conf_t conf = {0};
    bc_conf_init(&conf);

    static const struct option loptions[] = {
        {"mate"     , required_argument, NULL, 'm'},
        {"bc-start" , required_argument, NULL, 's'},
        {"bc-length", required_argument, NULL, 'l'},
        {"output"   , required_argument, NULL, 'o'},
        {"gz"       , required_argument, NULL, 'z'},
        {"help"     , no_argument      , NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    if (argc < 2) { usage(); return 1; }
    while ((c = getopt_long(argc, argv, ":l:m:o:s:h", loptions, NULL)) >= 0) {
        switch(c) {
            case 'l':
                conf.bc_length = (uint8_t)atoi(optarg);
                break;
            case 'm':
                conf.mate = (uint8_t)atoi(optarg);
                break;
            case 'o':
                conf.ofile = optarg;
                break;
            case 's':
                conf.bc_start = (int8_t)atoi(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case ':':
                usage();
                fprintf(stderr, "ERROR: option needs an argument: -%c\n", optopt);
                return 1;
            case '?':
                usage();
                fprintf(stderr, "ERROR: unrecognized option: -%c\n", optopt);
                return 1;
            default:
                usage();
                return 1;
        }
    }

    // Handle inputs
    if (conf.mate < 1 || conf.mate > 2) {
        usage();
        fprintf(stderr, "ERROR: -m,--mate must be 1 or 2\n");
        return 1;
    }

    if (conf.bc_start == 0) {
        fprintf(stderr, "ERROR: barcode start position should be 1-based, did you mean -s 1?\n");
        return 1;
    } else {
        conf.bc_start--;
    }

    // Init files and handle read errors
    if (optind >= argc) {
        usage();
        fprintf(stderr, "ERROR: no read FASTQ files provided\n");
        return 1;
    }

    gzFile fh1 = 0, fh2 = 0, oh1 = 0, oh2 = 0;
    kseq_t *ks1 = NULL, *ks2 = NULL;

    fh1 = gzopen(argv[optind], "r");
    if (!fh1) {
        fprintf(stderr, "ERROR: could not open input file: %s\n", argv[optind]);
        return 1;
    }
    ks1 = kseq_init(fh1);

    if (optind+1 < argc) {
        fh2 = gzopen(argv[optind+1], "r");
        if (!fh2) {
            fprintf(stderr, "ERROR: could not open input file: %s\n", argv[optind+1]);
            return 1;
        }
        ks2 = kseq_init(fh2);
    }

    if (conf.ofile) {
        oh1 = setup_output(conf.ofile, (ks2) ? 1 : 0);
        if (ks2) {
            oh2 = setup_output(conf.ofile, 2);
        }
    }

    // Process
    extract_barcodes(&conf, ks1, ks2, oh1, oh2);

    // Clean up
    if (oh2) { gzclose(oh2); }
    if (oh1) { gzclose(oh1); }
    if (ks2) { kseq_destroy(ks2); }
    if (ks1) { kseq_destroy(ks1); }
    if (fh2) { gzclose(fh2); }
    if (fh1) { gzclose(fh1); }

    return 0;
}
