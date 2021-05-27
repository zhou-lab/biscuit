/* generate QC files
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2021 Wanding.Zhou@vai.org
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

static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit qc [options] <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -h          This help\n");
    fprintf(stderr, "\n");
}

int main_qc(int argc, char *argv[]) {
    int c;

    if (argc < 1) { usage(); return 1; }
    while ((c = getopt(argc, argv, ":h")) >= 0) {
        switch (c) {
            case 'h': usage(); return 1;
            case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt);
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt);
            default: usage(); return 1;
        }
    }

    char *infn = optind < argc ? argv[optind++] : NULL;
    if (!infn) {
        usage();
        wzfatal("Please provide an input bam.\n");
    }

    fprintf(stdout, "You did it! You loaded %s\n", infn);

    return 0;
}
