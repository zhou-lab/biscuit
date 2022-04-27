/* generate QC files
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2021-2022 Jacob.Morrison@vai.org
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

void format_strand_report(FILE *fname, bsstrand_data_t *data) {
    fprintf(fname, "BISCUITqc Strand Table");
    int i;

    fprintf(fname, "\nStrand Distribution:\n");
    fprintf(fname, "strand\\BS      BSW (f)      BSC (r)\n");

    fprintf(fname, "     R1 (f):   ");
    for (i=0;i<2;++i) { fprintf(fname, "%-13d", data->strandcnt[i]); fprintf(fname, "\n"); }

    fprintf(fname, "     R1 (r):   ");
    for (i=0;i<2;++i) { fprintf(fname, "%-13d", data->strandcnt[4+i]); fprintf(fname, "\n"); }

    fprintf(fname, "     R2 (f):   ");
    for (i=0;i<2;++i) { fprintf(fname, "%-13d", data->strandcnt[8+i]); fprintf(fname, "\n"); }

    fprintf(fname, "     R2 (r):   ");
    for (i=0;i<2;++i) { fprintf(fname, "%-13d", data->strandcnt[12+i]); fprintf(fname, "\n"); }
}

void format_bsconv_report(FILE *fname, bsconv_data_t *data) {
    fprintf(fname, "BISCUITqc Conversion Rate by Read Average Table\n");
    fprintf(fname, "CpA\tCpC\tCpG\tCpT\n");

    int i;
    for (i=0; i<4; ++i) {
        if (i) { fprintf(fname, "\t"); }
        fprintf(fname, "%.8lf", (double) data->retn_conv_counts[2*i] / (data->retn_conv_counts[2*i] + data->retn_conv_counts[2*i+1]));
    }
    fprintf(fname, "\n");
}

void format_readpos_report(FILE *fname, cinread_data_t *data, char *type) {
    fprintf(fname, "BISCUITqc %s Retention by Read Position Table\n", type);
    fprintf(fname, "ReadInPair\tPosition\tConversion/Retention\tCount\n");

    unsigned i,j,k;
    char r;
    for (i=0; i<CIN_N_READS; ++i) {
        for (j=0; j<CIN_READ_LEN; ++j) {
            for (k=0; k<CIN_N_RET_STATES-1; ++k) { // ignore 'N' retention state
                if (k) {
                    r = 'R';
                } else {
                    r = 'C';
                }
                if (data->counts[i][j][k] > 0) {
                    fprintf(fname, "%u\t%u\t%c\t%d\n", i+1, j, r, data->counts[i][j][k]);
                }
            }
        }
    }
}

int process_qc(char *input_bam, qc_conf_t conf, bsstrand_data_t *data_bsstrand, bsconv_data_t *data_bsconv, cinread_data_t *data_cinread_cg, cinread_data_t *data_cinread_ch) {
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
        if (c->qual >= 40) {
            q40_tot++;
            cinread_func(b, 0, header, data_cinread_cg);
            cinread_func(b, 0, header, data_cinread_ch);
        }
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

            if (!(c->flag & BAM_FDUP) && (c->flag & BAM_FPAIRED) && (c->flag & BAM_FPROPER_PAIR) && (c->qual >= 40))
                bsconv_func(b, 0, header, data_bsconv);
        }

        // Determine bisulfite strand information
        bsstrand_func(b, 0, header, data_bsstrand);
    }
    bam_destroy1(b);

    sam_close(in);
    bam_hdr_destroy(header);

    if (ret != -1) /* truncated is -2 */
        wzfatal("[%s:%d] Alignment retrieval failed due to truncated file\n", __func__, __LINE__);

    format_mapq_table(conf.out.f_out_mapq, mapqs);
    format_dup_report(conf.out.f_out_dup, all_tot, all_dup, q40_tot, q40_dup);
    format_strand_report(conf.out.f_out_strand, data_bsstrand);
    format_bsconv_report(conf.out.f_out_bsconv, data_bsconv);
    format_readpos_report(conf.out.f_out_cgreadpos, data_cinread_cg, "CpG");
    format_readpos_report(conf.out.f_out_chreadpos, data_cinread_ch, "CpH");
    if (conf.single_end == 0)
        format_isize_table(conf.out.f_out_isize, isize, count_isizes);

    return ret;
}

static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit qc [options] <ref.fa> <in.bam> <sample_name>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -s    Run for single-end data\n");
    fprintf(stderr, "    -h    This help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note, this currently only produces a subset of QC metrics. Use scripts/QC.sh for full QC\n");
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
            case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt); break;
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt); break;
            default: usage(); return 1;
        }
    }

    char *reffn = optind < argc ? argv[optind++]: NULL;
    char *infn  = optind < argc ? argv[optind++] : NULL;
    char *samp  = optind < argc ? argv[optind++] : NULL;
    if (!reffn || !infn || !samp) {
        usage();
        wzfatal("Please provide a reference FASTA, input bam, or sample name.\n");
    }

    refcache_t *rs = init_refcache(reffn, 100, 100000);

    // bsstrand data and config
    bsstrand_conf_t conf_bsstrand = {0};
    bsstrand_data_t data_bsstrand = {0}; // set all counts to 0
    data_bsstrand.rs = rs;
    data_bsstrand.conf = &conf_bsstrand;

    // bsconv data and config
    bsconv_conf_t conf_bsconv = {0};
    conf_bsconv.max_cph = conf_bsconv.max_cpa = conf_bsconv.max_cpc = conf_bsconv.max_cpt = -1;
    conf_bsconv.max_cph_frac = 1.0; conf_bsconv.max_cpy_frac = 1.0;
    conf_bsconv.no_printing = 1;

    bsconv_data_t data_bsconv = {0};
    data_bsconv.rs = rs;
    data_bsconv.conf = &conf_bsconv;

    // CpG cinread data and config
    cinread_conf_t conf_cinread_cg = {0};
    conf_cinread_cg.skip_secondary = 1;
    conf_cinread_cg.skip_printing = 1;
    conf_cinread_cg.tgt = SL_CG;
    conf_cinread_cg.n_tp_names = 3;
    conf_cinread_cg.tp_names = malloc((conf_cinread_cg.n_tp_names)*sizeof(__tp_name_t));
    conf_cinread_cg.tp_names[0] = TP_QPAIR;
    conf_cinread_cg.tp_names[1] = TP_CQPOS;
    conf_cinread_cg.tp_names[2] = TP_CRETENTION;
    conf_cinread_cg.out = stdout;

    cinread_data_t data_cinread_cg = {0};
    data_cinread_cg.rs = rs;
    data_cinread_cg.conf = &conf_cinread_cg;

    // CpH cinread data and config
    cinread_conf_t conf_cinread_ch = {0};
    conf_cinread_ch.skip_secondary = 1;
    conf_cinread_ch.skip_printing = 1;
    conf_cinread_ch.tgt = SL_CH;
    conf_cinread_ch.n_tp_names = 3;
    conf_cinread_ch.tp_names = malloc((conf_cinread_ch.n_tp_names)*sizeof(__tp_name_t));
    conf_cinread_ch.tp_names[0] = TP_QPAIR;
    conf_cinread_ch.tp_names[1] = TP_CQPOS;
    conf_cinread_ch.tp_names[2] = TP_CRETENTION;
    conf_cinread_ch.out = stdout;

    cinread_data_t data_cinread_ch = {0};
    data_cinread_ch.rs = rs;
    data_cinread_ch.conf = &conf_cinread_ch;

    char out_mapq[LEN_SAMP+LEN_MAPQ+1];
    char out_dup[LEN_SAMP+LEN_DUP+1];
    char out_isize[LEN_SAMP+LEN_ISIZE+1];
    char out_strand[LEN_SAMP+LEN_STRAND+1];
    char out_bsconv[LEN_SAMP+LEN_BSCONV+1];
    char out_cgreadpos[LEN_SAMP+LEN_CGREADPOS+1];
    char out_chreadpos[LEN_SAMP+LEN_CHREADPOS+1];

    create_filename(out_mapq     , samp, "_mapq_table.txt");
    create_filename(out_dup      , samp, "_dup_report.txt");
    create_filename(out_strand   , samp, "_strand_table.txt");
    create_filename(out_bsconv   , samp, "_totalReadConversionRate.txt");
    create_filename(out_cgreadpos, samp, "_CpGRetentionByReadPos.txt");
    create_filename(out_chreadpos, samp, "_CpHRetentionByReadPos.txt");
    if (conf.single_end == 0)
        create_filename(out_isize, samp, "_isize_table.txt");
    
    qc_outfiles_t outfiles = {0};
    outfiles.f_out_mapq      = fopen(out_mapq  , "w");
    outfiles.f_out_dup       = fopen(out_dup   , "w");
    outfiles.f_out_strand    = fopen(out_strand, "w");
    outfiles.f_out_bsconv    = fopen(out_bsconv, "w");
    outfiles.f_out_cgreadpos = fopen(out_cgreadpos, "w");
    outfiles.f_out_chreadpos = fopen(out_chreadpos, "w");
    if (conf.single_end == 0)
        outfiles.f_out_isize = fopen(out_isize, "w");
    conf.out = outfiles;

    process_qc(infn, conf, &data_bsstrand, &data_bsconv, &data_cinread_cg, &data_cinread_ch);

    close_outfiles(outfiles, conf);
    free_refcache(rs);
    free(conf_cinread_cg.tp_names);
    free(conf_cinread_ch.tp_names);

    return 0;
}
