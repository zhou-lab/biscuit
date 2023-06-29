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

gzFile setup_output(const char *ofile) {
    // Length of file name, handle case where it's too short to do something with
    int len = strlen(ofile);
    if (len < 3) {
        fprintf(stderr, "ERROR: output file name not long enough to determine compression\n");
        return NULL;
    }

    // Default to writing with standard compression
    char *mode = "w6";

    // Change compression if not ".gz" file
    const char *ext = &ofile[len-3];
    if (strcmp(ext, ".gz") != 0) {
        mode = "w0";
    }
    fprintf(stdout, "mode: %s\n", mode);

    gzFile out = gzopen(ofile, mode);
    if (!out) {
        fprintf(stderr, "ERROR: could not open output file: %s\n", ofile);
        return NULL;
    }

    return out;
}

void extract_barcodes(kseq_t *ks1, kseq_t *ks2, gzFile oh1, gzFile oh2) {
    while (kseq_read(ks1) >= 0) {
        //fprintf(stdout, "%s", ks1->name.s);
        if (ks2 && kseq_read(ks2) < 0) { // read 2 has fewer reads
            fprintf(stderr, "WARNING: read 2 has fewer sequences\n");
            break;
        }

        //if (ks2) {
        //    fprintf(stdout, " %s", ks2->name.s);
        //}

        //fprintf(stdout, "\n");
    }
    if (ks2 && kseq_read(ks2) >= 0) {
        fprintf(stderr, "WARNING: read 1 has fewer sequences\n");
    }
}

static void usage() {
    bc_conf_t conf = {0};
    bc_conf_init(&conf);

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit bc [options] <FASTQ 1> [FASTQ 2]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Processing Options:\n");
    fprintf(stderr, "    -m, --mate INT         which mate the barcode is in (1 or 2) [%u]\n", conf.mate);
    fprintf(stderr, "    -s, --bc-start INT     start position of barcode in read (1-based) [%u]\n", conf.bc_start);
    fprintf(stderr, "    -l, --bc-length INT    length of barcode [%u]\n", conf.bc_length);
    fprintf(stderr, "Output Options:\n");
    fprintf(stderr, "    -o, --output1 STR      name of output file for read 1 (NULL writes to stdout) [NULL]\n");
    fprintf(stderr, "    -O, --output2 STR      name of output file for read 2 (NULL writes to stdout) [NULL]\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "    -h, --help             This help\n");
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
        {"output1"  , required_argument, NULL, 'o'},
        {"output2"  , required_argument, NULL, 'O'},
        {"help"     , no_argument      , NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    if (argc < 2) { usage(); return 1; }
    while ((c = getopt_long(argc, argv, ":l:m:o:s:O:h", loptions, NULL)) >= 0) {
        switch(c) {
            case 'l':
                conf.bc_length = (uint8_t)atoi(optarg);
                break;
            case 'm':
                conf.mate = (uint8_t)atoi(optarg);
                break;
            case 'o':
                conf.ofile1 = optarg;
                break;
            case 's':
                conf.bc_start = (int8_t)atoi(optarg);
                break;
            case 'O':
                conf.ofile2 = optarg;
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

    // Handle input errors
    if (conf.mate < 1 || conf.mate > 2) {
        usage();
        fprintf(stderr, "ERROR: -m,--mate must be 1 or 2\n");
        return 1;
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

    oh1 = setup_output(conf.ofile1);

    if (optind+1 < argc) {
        fh2 = gzopen(argv[optind+1], "r");
        if (!fh2) {
            fprintf(stderr, "ERROR: could not open input file: %s\n", argv[optind+1]);
            return 1;
        }
        ks2 = kseq_init(fh2);
        //oh2 = setup_output(conf.ofile2);
    }

    // Process
    extract_barcodes(ks1, ks2, oh1, oh2);

    // Clean up
    if (oh1) { gzclose(oh1); }
    if (ks2) { kseq_destroy(ks2); }
    if (ks1) { kseq_destroy(ks1); }
    //if (oh2) { gzclose(oh2); }
    //if (oh1) { gzclose(oh1); }
    if (fh2) { gzclose(fh2); }
    if (fh1) { gzclose(fh1); }

    return 0;
}
