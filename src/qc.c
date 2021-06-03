/* generate QC files
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2021 Jacob.Morrison@vai.org
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

#include "qc.h"

void format_mapq_table(FILE *fname, unsigned int *values) {
    fprintf(fname, "BISCUITqc Mapping Quality Table\n");
    fprintf(fname, "MapQ\tCount\n");

    fprintf(fname, "unmapped\t%u\n", values[N_MAPQ]);
    for (unsigned int i=0; i<N_MAPQ; i++)
        fprintf(fname, "%u\t%u\n", i, values[i]);
}

void format_isize_table(FILE *fname, unsigned int *values, unsigned int count) {
    fprintf(fname, "BISCUITqc Insert Size Table\n");
    fprintf(fname, "InsertSize\tFraction\tReadCount\n");

    for (unsigned int i=0; i<ISIZE+1; i++) {
        if (values[i] > 0) {
            double frac = values[i] / (double) count;
            fprintf(fname, "%u\t%.8lf\t%u\n", i, frac, values[i]);
        }
    }
}

void format_dup_report(FILE *fname, unsigned int all_tot, unsigned int all_dup, unsigned int q40_tot, unsigned int q40_dup) {
    fprintf(fname, "BISCUITqc Read Duplication Table\n");
    fprintf(fname, "Number of duplicate reads:\t%u\n", all_dup);
    fprintf(fname, "Number of reads:\t%u\n", all_tot);
    fprintf(fname, "Number of duplicate q40-reads:\t%u\n", q40_dup);
    fprintf(fname, "Number of q40-reads:\t%u\n", q40_tot);
}

int process_qc(char *input_bam, qc_conf_t conf) {
    int ret = 0;
    samFile *in = sam_open(input_bam, "rb");
    bam_hdr_t *header = sam_hdr_read(in);

    unsigned int all_tot  = 0; // Number of reads in file
    unsigned int all_dup  = 0; // Number of duplicate flagged reads
    unsigned int q40_tot  = 0; // Number of reads in file with MAPQ >= 40
    unsigned int q40_dup  = 0; // Number of duplicate flagged reads with MAPQ >= 40
    unsigned int count_isizes    = 0; // Number of reads with insert size
    unsigned int mapqs[N_MAPQ+1] = {0}; // Stores MAPQ values
    unsigned int isize[ISIZE+1]  = {0}; // Stores insert size values
    bam1_t *b = bam_init1();
    while ((ret = sam_read1(in, header, b)) >= 0) {
        bam1_core_t *c = &b->core;

        // Counts values for duplicate rates
        all_tot++;
        if (c->flag & BAM_FDUP)
            all_dup++;
        if (c->qual >= 40)
            q40_tot++;
        if ((c->flag & BAM_FDUP) && (c->qual >= 40))
            q40_dup++;

        // Retrieve MAPQ and insert size values
        if (!(c->flag & BAM_FSECONDARY)) {
            if (c->flag & BAM_FUNMAP)
                mapqs[N_MAPQ]++;
            else
                mapqs[c->qual]++;
        
            if ((conf.single_end == 0) && (c->flag & BAM_FPROPER_PAIR) && (c->qual >= 40)) {
                if ((c->isize >= 0) && (c->isize <= ISIZE)) {
                    count_isizes++;
                    isize[c->isize]++;
                }
            }
        }
    }
    bam_destroy1(b);

    sam_close(in);
    bam_hdr_destroy(header);

    if (ret != -1) /* truncated is -2 */
        wzfatal("[%s:%d] Alignment retrieval failed due to truncated file\n", __func__, __LINE__);

    format_mapq_table(conf.out.f_out_mapq, mapqs);
    format_dup_report(conf.out.f_out_dup, all_tot, all_dup, q40_tot, q40_dup);
    if (conf.single_end == 0)
        format_isize_table(conf.out.f_out_isize, isize, count_isizes);

    return ret;
}

// TODO: May need to change the name as this moves beyond just alignment QC
static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit qc [options] <in.bam> <sample_name>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -s    Run for single-end data\n");
    fprintf(stderr, "    -h    This help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note, this currently only produces alignment QC metrics. Use scripts/QC.sh for full QC\n");
    fprintf(stderr, "\n");
}

int main_qc(int argc, char *argv[]) {
    int c;
    qc_conf_t conf = {0};
    conf.single_end = 0;

    if (argc < 2) { usage(); return 1; }
    while ((c = getopt(argc, argv, ":hs")) >= 0) {
        switch (c) {
            case 's': conf.single_end = 1; break;
            case 'h': usage(); return 1;
            case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt);
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt);
            default: usage(); return 1;
        }
    }

    char *infn = optind < argc ? argv[optind++] : NULL;
    char *samp = optind < argc ? argv[optind++] : NULL;
    if (!infn || !samp) {
        usage();
        wzfatal("Please provide an input bam or sample name.\n");
    }

    char out_mapq[LEN_SAMP+LEN_DUP+1];
    char out_dup[LEN_SAMP+LEN_DUP+1];
    char out_isize[LEN_SAMP+LEN_ISIZE+1];

    create_filename(out_mapq , samp, "_mapq_table.txt");
    create_filename(out_dup  , samp, "_dup_report.txt");
    if (conf.single_end == 0)
        create_filename(out_isize, samp, "_isize_table.txt");
    
    qc_outfiles_t outfiles = {0};
    outfiles.f_out_mapq  = fopen(out_mapq , "w");
    outfiles.f_out_dup   = fopen(out_dup  , "w");
    if (conf.single_end == 0)
        outfiles.f_out_isize = fopen(out_isize, "w");
    conf.out = outfiles;

    fprintf(stdout, "You did it! You loaded %s\n", infn);
    process_qc(infn, conf);

    close_outfiles(outfiles, conf);

    return 0;
}
