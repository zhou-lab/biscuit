/* convert bam to epiread format with supplied snp bed file
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
 */
#include "epiread.h"

const char SKIP_EPI = '-';
const char SKIP_INS = 'i';
const char SKIP_DEL = 'd';
const char FILTERED = 'F';
const char IGNORED  = 'x';
const char DELETION = 'D';
const char SOFTCLIP = 'P';
const char METHYLAT = 'M';
const char UNMETHYL = 'U';
const char OPEN_ACC = 'O';
const char SHUT_ACC = 'S';
const char AMBIG_GA = 'R';
const char AMBIG_CT = 'Y';

// TODO: Work these variables into two variables in the conf_t type and can be set via CLI options
// For short reads
#define MAX_READ_LENGTH 302
#define MAX_RLEN 50

// For long reads, I'm guessing there's a better way to do this, but until I
// know what that is, this is going to be the way it's done
// This also probably too long for PacBio and could be too short for ONT, but
// I'll put in some checks to skip lines that are too long
#define MAX_READ_LENGTH_LR 50000
#define MAX_RLEN_LR 500

DEFINE_VECTOR(int_v, int)
DEFINE_VECTOR(char_v, char)

typedef struct episnp_chrom1_t {
    char     *chrm; /* chromosome */
    size_t    n;    /* number of snps */
    uint32_t *locs; /* snp locations */
    uint8_t  *meth; /* is methylation callable for location */
} episnp_chrom1_t;

DEFINE_VECTOR(episnp_chrom1_v, episnp_chrom1_t)

void destroy_episnp(episnp_chrom1_v *episnp) {
    uint32_t i;
    for (i=0; i<episnp->size; ++i) {
        episnp_chrom1_t *e = ref_episnp_chrom1_v(episnp, i);
        free(e->locs);
        free(e->meth);
        free(e->chrm);
    }
    free_episnp_chrom1_v(episnp);
}

// Get all episnps from one chromosome
static inline episnp_chrom1_t *get_episnp1(episnp_chrom1_v *episnp, char *chrm) {
    uint32_t i;
    episnp_chrom1_t *episnp1;
    for (i=0; i<episnp->size; ++i) {
        episnp1 = ref_episnp_chrom1_v(episnp, i);
        if (strcmp(episnp1->chrm, chrm) == 0) return episnp1;
    }
    return NULL;
}

static inline episnp_chrom1_t *get_n_insert_episnp1(episnp_chrom1_v *episnp, char *chrm) {
    episnp_chrom1_t *episnp1 = get_episnp1(episnp, chrm);
    if (!episnp1) {
        episnp1 = next_ref_episnp_chrom1_v(episnp);
        episnp1->chrm = strdup(chrm);
        episnp1->locs = NULL;
        episnp1->meth = NULL;
        episnp1->n = 0;
    }
    return episnp1;
}

#define episnp_test(snps, i) snps[(i)>>3]&(1<<((i)&0x7))
#define episnp_set(snps, i) snps[(i)>>3] |= 1<<((i)&0x7)

typedef struct {
    char             *bam_fn; /* BAM filename */
    char             *ref_fn; /* reference filename */
    episnp_chrom1_v  *snp;    /* vector of snp locations */
    wqueue_t(window) *q;      /* queue window */
    wqueue_t(record) *rq;     /* queue records */
    epiread_conf_t           *conf;   /* config variables */
} result_t;

typedef struct {
    wqueue_t(record) *q;
    int n_bams;
    char **bam_fns;
    char *outfn;
    char *statsfn;
    char *header;
    target_v *targets;
    epiread_conf_t *conf;
} writer_conf_t;

static void *epiread_write_func(void *data) {
    writer_conf_t *c = (writer_conf_t*) data;

    FILE *out;
    if (c->outfn) {
        out = fopen(c->outfn, "w");
    } else {
        out = stdout;
    }

    int64_t next_block = 0;
    record_v *records = init_record_v(20);

    while (1) {
        record_t rec;
        wqueue_get(record, c->q, &rec);
        if(rec.block_id == RECORD_QUEUE_END) break;

        if (rec.block_id == next_block) {
            do {
                if (rec.s.s) fputs(rec.s.s, out);
                free(rec.s.s);

                // Get next block from shelf if available else return OBSOLETE and retrieve new block from queue
                next_block++;
                pop_record_by_block_id(records, next_block, &rec);
            } while (rec.block_id != RECORD_SLOT_OBSOLETE);
        } else {
            // Shelf the block if not next
            put_into_record_v(records, rec);
        }
    }

    free_record_v(records);
    if (c->outfn) {
        // For stdout, will close at the end of main
        fflush(out);
        fclose(out);
    }

    return 0;
}

void run_length_encode(char *str, char *out, epiread_conf_t *conf) {
    int run_len;
    char *count;
    if (conf->is_long_read) {
        count = (char *)calloc(MAX_RLEN_LR, sizeof(char));
    } else {
        count = (char *)calloc(MAX_RLEN, sizeof(char));
    }
    int len = strlen(str);

    int i, j=0, k;

    for (i=0; i<len; i++) {
        out[j++] = str[i];

        run_len = 1;
        while (i+1 < len && str[i] == str[i+1]) {
            run_len++;
            i++;
        }

        // Currently, the code is set up to set ABC = ABC. If you want ABC = A1B1C1, then comment out the if-statement,
        // but leave the sprintf and for loop untouched.
        if (run_len > 1) {
            sprintf(count, "%d", run_len);
            for (k=0; *(count+k); k++, j++) {
                out[j] = count[k];
            }
        }
    }

    free(count);

    out[j] = '\0';
}

// Format one bam record into the epi-bed format (positions are 0-based)
// TODO: This function is getting verbose, refactor/cleanup?
static void format_epi_bed(
        kstring_t *epi, bam1_t *b, uint8_t bsstrand, char *chrm, window_t *w, epiread_conf_t *conf,
        char *rle_arr_cg, char *rle_arr_gc, char *rle_arr_vr, uint32_t start, uint32_t end) {

    // Set max read length
    uint32_t max_read_length = (conf->is_long_read) ? MAX_READ_LENGTH_LR : MAX_READ_LENGTH;

    // Only shift the beginning of the window if conf->epiread_reg_start == w->beg (only occurs for the first window)
    // This catches reads that start before the first window when calling epiread with a region (-g option) provided
    uint32_t print_w_beg = conf->epiread_reg_start == w->beg ? w->beg - max_read_length : w->beg;

    // Only shift the end of the window if conf->epiread_reg_end == w->end (only occurs for the last window)
    // This catches reads that start at the last location when calling epiread with a region (-g option) provided
    uint32_t print_w_end = conf->epiread_reg_end == w->end ? w->end + max_read_length : w->end;

    // Columns: chromosome, start, end, read name, read number, BS strand, encoded CG RLE
    // If running in NOMe-seq mode, then encoded GC RLE is added as a last column
    if (start > 0 && (unsigned) start >= print_w_beg && (unsigned) start < print_w_end) {
        uint8_t write_read_cg = 1;
        uint8_t write_read_gc = 1;
        uint8_t write_read_vr = 1;
        int len_cg = (int)strlen(rle_arr_cg);
        int len_vr = (int)strlen(rle_arr_vr);
        int len_gc = conf->comm.is_nome ? (int)strlen(rle_arr_gc) : -1;

        // Check if RLE string is only composed of F's and x's and P's
        if (conf->filter_empty_epiread) {
            char *filt_ignr_soft = "FxP";

            size_t check_rle_cg = strspn(rle_arr_cg, filt_ignr_soft);
            if (len_cg == (int)check_rle_cg) { write_read_cg = 0; }

            size_t check_rle_vr = strspn(rle_arr_vr, filt_ignr_soft);
            if (len_vr == (int)check_rle_vr) { write_read_vr = 0; }

            if (conf->comm.is_nome) {
                size_t check_rle_gc = strspn(rle_arr_gc, filt_ignr_soft);
                if (len_gc == (int)check_rle_gc) { write_read_gc = 0; }
            } else {
                write_read_gc = 0;
            }
        }

        if (write_read_cg || write_read_gc || write_read_vr) {
            ksprintf(epi, "%s\t%d\t%d\t%s\t%c\t%c",
                    chrm,
                    start-1,
                    end,
                    bam_get_qname(b),
                    (b->core.flag&BAM_FREAD2) ? '2' : '1',
                    bsstrand ? '-' : '+');

            // Encoded string can be no larger than 2 times the original string length
            char *encoded_rle = (char *) malloc(sizeof(char) * (2*len_cg+1));
            run_length_encode(rle_arr_cg, encoded_rle, conf);
            ksprintf(epi, "\t%s", encoded_rle);

            if (conf->comm.is_nome) {
                // Encoded string can be no larger than 2 times the original string length
                char *encoded_rle_gc = (char *) malloc(sizeof(char) * (2*len_gc+1));
                run_length_encode(rle_arr_gc, encoded_rle_gc, conf);
                ksprintf(epi, "\t%s", encoded_rle_gc);
                free(encoded_rle_gc);
            } else {
                ksprintf(epi, "\t.");
            }

            char *encoded_rle_vr = (char *) malloc(sizeof(char) * (2*len_vr+1));
            run_length_encode(rle_arr_vr, encoded_rle_vr, conf);
            ksprintf(epi, "\t%s", encoded_rle_vr);

            // end line
            kputc('\n', epi);

            free(encoded_rle_vr);
            free(encoded_rle);
        } else {
            if (conf->comm.verbose) {
                fprintf(stderr, "Filtering CG read: %s\n", rle_arr_cg);
                fprintf(stderr, "Filtering GC read: %s\n", rle_arr_gc);
                fprintf(stderr, "Filtering variant read: %s\n", rle_arr_vr);
            }
        }
    }
}

// Format one bam record into the old epiread format (positions are 0-based)
// TODO: This function is getting verbose, refactor/cleanup?
static void format_epiread_old(
        kstring_t *epi, bam1_t *b, uint8_t bsstrand, char *chrm, window_t *w, uint8_t *snps, epiread_conf_t *conf,
        int_v *snp_p, int_v *hcg_p, int_v *gch_p, int_v *cg_p,
        char_v *snp_c, char_v *hcg_c, char_v *gch_c, char_v *cg_c) {

    // Set max read length
    uint32_t max_read_length = (conf->is_long_read) ? MAX_READ_LENGTH_LR : MAX_READ_LENGTH;

    // Only shift the beginning of the window if conf->epiread_reg_start == w->beg (only occurs for the first window)
    // This catches reads that start before the first window when calling epiread with a region (-g option) provided
    uint32_t print_w_beg = conf->epiread_reg_start == w->beg ? w->beg - max_read_length : w->beg;

    // Only shift the end of the window if conf->epiread_reg_end == w->end (only occurs for the last window)
    // This catches reads that start at the last location when calling epiread with a region (-g option) provided
    uint32_t print_w_end = conf->epiread_reg_end == w->end ? w->end + max_read_length : w->end;

    uint32_t k;

    if (conf->comm.is_nome) { // nome-seq
        // TODO: This solution has the same issues as those in the bsseq option in the else-block of this
        //       if-statement
        int first_epi = 0;
        if (hcg_p->size > 0 && gch_p->size > 0) {
            first_epi = get_int_v(hcg_p, 0) < get_int_v(gch_p, 0) ? get_int_v(hcg_p, 0) : get_int_v(gch_p, 0);
        } else if (hcg_p->size > 0 && gch_p->size == 0) {
            first_epi = get_int_v(hcg_p, 0);
        } else if (hcg_p->size == 0 && gch_p->size > 0) {
            first_epi = get_int_v(gch_p, 0);
        }

        // Avoid double counting between windows
        if (first_epi > 0 && (unsigned) first_epi >= print_w_beg && (unsigned) first_epi < print_w_end) {
            ksprintf(epi, "%s\t%s\t%c\t%c",
                    chrm,
                    bam_get_qname(b),
                    (b->core.flag&BAM_FREAD2) ? '2' : '1',
                    bsstrand ? '-' : '+');

            // HCG context (0-based)
            if (hcg_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(hcg_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<hcg_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(hcg_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(hcg_c, 0));
                for (k=1; k<hcg_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(hcg_c, k));
            } else {
                kputs("\t.\t.", epi);
            }

            // GCH context (0-based)
            if (gch_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(gch_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<gch_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(gch_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(gch_c, 0));
                for (k=1; k<gch_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(gch_c, k));
            } else {
                kputs("\t.\t.", epi);
            }

            // snp (0-based)
            if (snp_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(snp_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<snp_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(snp_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(snp_c, 0));
                for (k=1; k<snp_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(snp_c, k));
            } else if (snps) {
                kputs("\t.\t.", epi);
            } else {
                kputs("\t\t", epi);
            }

            // end line
            kputc('\n', epi);
        }
    } else { // bs-seq
        // Avoid double counting between windows
        // TODO: While this removes the problem with depending on a potentially unitialized value in the
        //       if-statement, it means that any time we have a cg vector that is empty, we won't print
        //       anything. I think the desired output is that empty CGs will have "\t.\t." printed, which
        //       I don't think will happen in this fix. As this is the old format, I'm not too concerned
        //       about getting the same output as this would ideally be deprecated/removed in the future.
        int cg_start = cg_p->size > 0 ? get_int_v(cg_p, 0) : 0;
        if (cg_start > 0 && (unsigned) cg_start >= print_w_beg && (unsigned) cg_start < print_w_end) {
            ksprintf(epi, "%s\t%s\t%c\t%c",
                    chrm,
                    bam_get_qname(b),
                    (b->core.flag&BAM_FREAD2) ? '2' : '1',
                    bsstrand ? '-' : '+');

            // CpG context (0-based)
            if (cg_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(cg_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<cg_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(cg_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(cg_c, 0));
                for (k=1; k<cg_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(cg_c, k));
            } else {
                kputs("\t.\t.", epi);
            }

            // snp (0-based)
            if (snp_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(snp_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<snp_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(snp_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(snp_c, 0));
                for (k=1; k<snp_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(snp_c, k));
            } else if (snps) {
                kputs("\t.\t.", epi);
            } else {
                kputs("\t\t", epi);
            }

            // end line
            kputc('\n', epi);
        }
    }
}

// Format one bam record to pairwise format (positions are 1-based)
// TODO: This function is getting verbose, refactor/cleanup?
static void format_epiread_pairwise(
        kstring_t *epi, char *chrm, window_t *w, epiread_conf_t *conf,
        int_v *snp_p, int_v *hcg_p, int_v *gch_p, int_v *cg_p,
        char_v *snp_c, char_v *hcg_c, char_v *gch_c, char_v *cg_c) {

    // Set max read length
    uint32_t max_read_length = (conf->is_long_read) ? MAX_READ_LENGTH_LR : MAX_READ_LENGTH;

    // Only shift the beginning of the window if conf->epiread_reg_start == w->beg (only occurs for the first window)
    // This catches reads that start before the first window when calling epiread with a region (-g option) provided
    uint32_t print_w_beg = conf->epiread_reg_start == w->beg ? w->beg - max_read_length : w->beg;

    // Only shift the end of the window if conf->epiread_reg_end == w->end (only occurs for the last window)
    // This catches reads that start at the last location when calling epiread with a region (-g option) provided
    uint32_t print_w_end = conf->epiread_reg_end == w->end ? w->end + max_read_length : w->end;

    uint32_t j, k;
    for (k=0; k<snp_p->size; ++k) {
        // avoid double counting between windows
        if (!((unsigned) get_int_v(snp_p, k) >= print_w_beg && (unsigned) get_int_v(snp_p, k) < print_w_end))
            continue;

        if (conf->comm.is_nome) { // nome-seq
            for (j=0; j<hcg_p->size; ++j) { // snp and HCG context
                ksprintf(epi, "%s\t%d\t%d\t%c\t%c\n",
                        chrm,
                        get_int_v(snp_p, k),
                        get_int_v(hcg_p, j),
                        get_char_v(snp_c, k),
                        get_char_v(hcg_c, j));
            }
            for (j=0; j<gch_p->size; ++j) { // snp and GCH context
                ksprintf(epi, "%s\t%d\t%d\t%c\t%c\n",
                        chrm,
                        get_int_v(snp_p, k),
                        get_int_v(gch_p, j),
                        get_char_v(snp_c, k),
                        get_char_v(gch_c, j));
            }
        } else { // bs-seq
            for (j=0; j<cg_p->size; ++j) {
                // chrm, snp position, cpg position, snp calling, cytosine calling
                ksprintf(epi, "%s\t%d\t%d\t%c\t%c\n",
                        chrm,
                        get_int_v(snp_p, k),
                        get_int_v(cg_p, j),
                        get_char_v(snp_c, k),
                        get_char_v(cg_c, j));
            }
        }
    }
}

void skipped_base_old(
        refcache_t *rs, char rb, uint8_t bss, uint32_t rj, uint32_t qj, epiread_conf_t *conf, char skip_epi,
        int_v *hcg_p, int_v *gch_p, int_v *cg_p, char_v *hcg_c, char_v *gch_c, char_v *cg_c) {
    if (bss && rb == 'G' && rj-1 >= rs->beg) {
        char rb0 = refcache_getbase_upcase(rs, rj-1);
        if (conf->comm.is_nome) {
            if (rj+1 <= rs->end) {
                char rb1 = refcache_getbase_upcase(rs, rj+1);
                if (rb0 == 'C' && rb1 != 'C' && qj > 0) {
                    push_int_v(hcg_p, (int) rj-1); push_char_v(hcg_c, skip_epi);
                } else if (rb0 != 'C' && rb1 == 'C') {
                    push_int_v(gch_p, (int) rj); push_char_v(gch_c, skip_epi);
                }
            }
        } else {
            if (rb0 == 'C') {
                push_int_v(cg_p, (int) rj-1); push_char_v(cg_c, skip_epi);
            }
        }
    }
    if (!bss && rb == 'C' && rj+1 <= rs->end) {
        char rb1 = refcache_getbase_upcase(rs, rj+1);
        if (conf->comm.is_nome) {
            if (rj-1 >= rs->beg) {
                char rb0 = refcache_getbase_upcase(rs, rj-1);
                if (rb0 != 'G' && rb1 == 'G') {
                    push_int_v(hcg_p, (int) rj); push_char_v(hcg_c, skip_epi);
                } else if (rb0 == 'G' && rb1 != 'G') {
                    push_int_v(gch_p, (int) rj); push_char_v(gch_c, skip_epi);
                }
            }
        } else {
            if (rb1 == 'G') {
                push_int_v(cg_p, (int) rj); push_char_v(cg_c, skip_epi);
            }
        }
    }
}

static inline void add_filtered(char *cg, char *var, char *gc, uint32_t idx) {
    cg[idx]  = FILTERED;
    var[idx] = FILTERED;
    gc[idx]  = FILTERED;
}

static void *process_func(void *data) {
    result_t *res  = (result_t*) data;
    epiread_conf_t   *conf = (epiread_conf_t*) res->conf;

    htsFile   *in  = hts_open(res->bam_fn, "rb");
    hts_idx_t *idx = sam_index_load(in, res->bam_fn);
    if (!idx) {
        fprintf(stderr, "[%s:%d] BAM %s is not indexed?\n", __func__, __LINE__, res->bam_fn);
        fflush(stderr);
        exit(1);
    }
    bam_hdr_t *header = sam_hdr_read(in);

    refcache_t *rs = init_refcache(res->ref_fn, 1000, 1000);
    uint32_t j;

    record_t rec;
    memset(&rec, 0, sizeof(record_t));
    window_t w;
    while (1) {
        wqueue_get(window, res->q, &w);
        if (w.tid == -1) break;

        rec.tid = w.tid;
        char *chrm = header->target_name[w.tid];

        uint32_t snp_beg = w.beg>1000 ? w.beg-1000 : 1; // start location of snps
        uint32_t snp_end = w.end+1000;

        // Make snp lookup table (if supplied)
        uint8_t *snps = NULL;
        uint8_t *meth = NULL;
        if (res->snp) {
            snps = calloc((snp_end-snp_beg)/8+1, sizeof(uint8_t));
            meth = calloc((snp_end-snp_beg)/8+1, sizeof(uint8_t));
            episnp_chrom1_t *episnp1 = get_episnp1(res->snp, chrm);

            // If chromosome is found in snp file
            if (episnp1) {
                for (j=0; j<episnp1->n; ++j) {
                    uint32_t l = episnp1->locs[j];
                    uint8_t  m = episnp1->meth[j];
                    if (l>=snp_beg && l<snp_end) {
                        episnp_set(snps, l-snp_beg);
                        if (m) { episnp_set(meth, l-snp_beg); }
                    }
                }
            }
        }

        // The epiread string
        rec.s.l = rec.s.m = 0; rec.s.s = 0;

        refcache_fetch(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
        hts_itr_t *iter = sam_itr_queryi(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
        bam1_t *b = bam_init1();
        hts_base_mod_state *mod_state = NULL;
        hts_base_mod mod[5] = {0};
        int cutoff_score = (int)(conf->modbam_prob * 255);
        int mod_pos, n_mods;
        if (conf->use_modbam) {
            mod_state = hts_base_mod_state_alloc();
        }
        int ret;

        while ((ret = sam_itr_next(in, iter, b))>0) {
            bam1_core_t *c = &b->core;

            // Set up modBAM tags extraction
            if (conf->use_modbam) {
                if (bam_parse_basemod2(b, mod_state, HTS_MOD_REPORT_UNCHECKED) < 0) {
                    wzfatal("ERROR: Failed to parse base modifications\n");
                }

                // Check for correct number of modifications
                int n_all_mods = 0;
                int *all_mods = bam_mods_recorded(mod_state, &n_all_mods);
                if (n_all_mods > 1) {
                    wzfatal("ERROR: too many modifications found. Only one modification allowed per read.\n");
                }

                // Already ensured we only have one modification (index = 0), check for correct mod type now
                int m_strand, m_implicit;
                char m_canonical;
                bam_mods_queryi(mod_state, 0, &m_strand, &m_implicit, &m_canonical);
                if (all_mods[0] != 'm') {
                    wzfatal("ERROR: must be a methylation modification ('m')\n");
                }
                if (m_canonical != 'C' && m_canonical != 'G') {
                    wzfatal("ERROR: modification must fall on a C or G\n");
                }
            }

            // Read-based filtering
            if (c->qual < conf->filt.min_mapq) continue;
            if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->filt.min_read_len) continue;
            if (c->flag > 0) { // only when any flag is set
                if (conf->filt.filter_secondary && c->flag & BAM_FSECONDARY) continue;
                if (conf->filt.filter_duplicate && c->flag & BAM_FDUP) continue;
                if (conf->filt.filter_ppair && c->flag & BAM_FPAIRED && !(c->flag & BAM_FPROPER_PAIR)) continue;
                if (conf->filt.filter_qcfail && c->flag & BAM_FQCFAIL) continue;
            }

            uint8_t *nm = bam_aux_get(b, "NM");
            if (nm && bam_aux2i(nm) > conf->filt.max_nm) continue;

            uint8_t *as = bam_aux_get(b, "AS");
            if (as && bam_aux2i(as) < conf->filt.min_score) continue;

            // For now, assume if using modBAM tags that you're running ONT data that has no concept
            // of a bisulfite strand or retained cytosines
            uint8_t bsstrand = conf->use_modbam ? 0 : get_bsstrand(rs, b, conf->filt.min_base_qual, 0);
            uint32_t cnt_ret = conf->use_modbam ? 0 : cnt_retention(rs, b, bsstrand);
            if (cnt_ret > conf->filt.max_retention) continue;

            // Pairwise epiread format variables
            int_v  *snp_p = init_int_v(10);  // snp position
            char_v *snp_c = init_char_v(10); // snp character
            int_v  *cg_p=0, *hcg_p=0, *gch_p=0;
            char_v *cg_c=0, *hcg_c=0, *gch_c=0;
            if (conf->comm.is_nome) {
                hcg_p = init_int_v(10);  // hcg positions
                hcg_c = init_char_v(10); // hcg characters
                gch_p = init_int_v(10);  // gch positions
                gch_c = init_char_v(10); // gch characters
            } else {
                cg_p = init_int_v(10);  // cpg positions
                cg_c = init_char_v(10); // cpg characters
            }

            // Run length encoding strings
            char *rle_arr_cg; // use qpos+j to determine which index will be written
            char *rle_arr_gc;
            char *rle_arr_vr;
            if (conf->is_long_read) {
                rle_arr_cg = (char *)calloc(MAX_READ_LENGTH_LR, sizeof(char));
                rle_arr_vr = (char *)calloc(MAX_READ_LENGTH_LR, sizeof(char));
                rle_arr_gc = (char *)calloc(MAX_READ_LENGTH_LR, sizeof(char));
            } else {
                rle_arr_cg = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
                rle_arr_vr = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
                rle_arr_gc = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
            }

            uint32_t i, j;
            uint8_t rle_set        = 0; // Know when to write info to array generally
            uint32_t n_deletions    = 0; // Number of deletions, used to shift the RLE string accordingly when deletions are present
            uint32_t n_insertions   = 0; // Number of insertions, used to shift the end position of the RLE string when insertions are present
            uint32_t softclip_start = 0; // Number of soft clip bases occurring at start of read - used to adjust starting position

            uint32_t rpos  = c->pos  + 1; // 1-based reference position
            uint32_t rmpos = c->mpos + 1; // 1-based mate reference position
            uint32_t qpos  = 0;           // query position
            uint32_t qj    = 0;           // current position in query
            uint32_t qjd   = 0;           // current position in query adjusted by n_deletions

            // Handle read lengths for overlapping reads
            uint32_t read_length = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));

            uint32_t mate_length;
            uint8_t *mc = bam_aux_get(b, "MC");
            if (mc) {
                mc++;
                mate_length = get_mate_length((char *)mc);
            } else {
                // If MC tag is missing, then assume reads are the same length
                mate_length = read_length;
            }

            // -1 accounts for the 1-based nature of the coordinates
            uint32_t rend  = rpos  + read_length - 1;
            uint32_t rmend = rmpos + mate_length - 1;

            char rb, qb;
            for (i=0; i<c->n_cigar; ++i) {
                uint32_t op    = bam_cigar_op(bam_get_cigar(b)[i]);
                uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                switch(op) {
                    case BAM_CMATCH:
                    case BAM_CEQUAL:
                    case BAM_CDIFF:
                        for (j=0; j<oplen; ++j) {
                            // Query positions
                            qj = qpos + j;
                            qjd = qj + n_deletions;

                            // Reference and query bases
                            rb = refcache_getbase_upcase(rs, rpos+j);
                            qb = bscall(b, qj);

                            // Get mod tags before filtering to make sure we properly count deltas in MM tag
                            if (conf->use_modbam) {
                                n_mods = bam_mods_at_next_pos(b, mod_state, mod, 1);
                                if (n_mods < 0) {
                                    wzfatal("ERROR: problem encountered retrieving next base modification\n");
                                }
                            }

                            // Base filtering
                            // Low quality bases
                            if (bam_get_qual(b)[qj] < conf->filt.min_base_qual) {
                                skipped_base_old(rs, rb, bsstrand, rpos+j, qj, conf, SKIP_EPI, hcg_p, gch_p, cg_p, hcg_c, gch_c, cg_c);
                                add_filtered(rle_arr_cg, rle_arr_vr, rle_arr_gc, qjd);

                                continue;
                            }

                            // Read-position-based filtering
                            // Follows the same form as pileup, so qpos is adjusted to be 1-based/1-indexed and
                            // filtering is done according to that, rather than being 0-based/0-indexed
                            if (qj+1 <= conf->filt.min_dist_end_5p || c->l_qseq < (int32_t)(qj+1 + conf->filt.min_dist_end_3p)) {
                                skipped_base_old(rs, rb, bsstrand, rpos+j, qj, conf, SKIP_EPI, hcg_p, gch_p, cg_p, hcg_c, gch_c, cg_c);
                                add_filtered(rle_arr_cg, rle_arr_vr, rle_arr_gc, qjd);

                                continue;
                            }

                            /* If reads 1 and 2 overlap, skip counting bases in read 2
                             * Read lengths are relative to the reference as defined by the CIGAR string (if MC tag
                             *     not given in read, then mate length is assumed to be the same as the current
                             *     read)
                             * Overlapping bases in both proper and improper pairs (if user requests these be
                             *     included) will be ignored
                             *
                             * The filtering removes bases from read 2 (usually the read on the complement strand)
                             * that fall into the overlapped region.
                             */
                            if ((conf->filt.filter_doublecnt) && (c->flag & BAM_FREAD2) &&
                                (rpos+j >= max(rpos, rmpos)) && (rpos+j <= min(rend, rmend))) {
                                skipped_base_old(rs, rb, bsstrand, rpos+j, qj, conf, SKIP_EPI, hcg_p, gch_p, cg_p, hcg_c, gch_c, cg_c);
                                add_filtered(rle_arr_cg, rle_arr_vr, rle_arr_gc, qjd);

                                continue;
                            }

                            // Methylation is handled differently for modBAMs and regular BAMs
                            if (conf->use_modbam) {
                                uint8_t is_cpg = is_modbam_cpg(c->flag, mod[0].strand, mod[0].canonical_base, qb, rb, rs, rpos+j);

                                if (conf->use_modbam && n_mods > 0) {
                                    float mod_probability = calculate_mod_probability(mod[0].qual);
                                    //fprintf(stderr, "pos: %i, qj: %i, can base: %c, qb: %c, n_mods: %i, qual: %i, prob: %f, is_cpg: %u\n",
                                    //        c->pos+1 - softclip_start + qjd - n_insertions,
                                    //        qj, mod[0].canonical_base, qb, n_mods, mod[0].qual, mod_probability, is_cpg);
                                    push_int_v(cg_p, (int) rpos+j);
                                    if (is_cpg && mod[0].qual >= 0 && mod_probability > conf->modbam_prob) {
                                        push_char_v(cg_c, 'C');
                                        rle_arr_cg[qjd] = METHYLAT;
                                        rle_set = 1;
                                    } else if (is_cpg && mod[0].qual >= 0 && mod_probability < 1.0-conf->modbam_prob) {
                                        push_char_v(cg_c, 'T');
                                        rle_arr_cg[qjd] = UNMETHYL;
                                        rle_set = 1;
                                    } else {
                                        push_char_v(cg_c, 'N');
                                    }
                                }
                            } else {
                                // reference is a G
                                if (bsstrand && rb == 'G' && rpos+j-1 >= rs->beg) {
                                    char rb0 = refcache_getbase_upcase(rs, rpos+j-1); // previous base
                                    if (conf->comm.is_nome) { // nome-seq
                                        if (rpos+j+1 <= rs->end) { // prevent overflow
                                            char rb1 = refcache_getbase_upcase(rs, rpos+j+1); // next base
                                            if (rb0 == 'C' && rb1 != 'C') { // HCG context
                                                // Note: measure G in CpG context, record location of C
                                                if (qj > 0) { push_int_v(hcg_p, (int) rpos+j-1); }
                                                if (qb == 'A') {
                                                    push_char_v(hcg_c, 'T');
                                                    rle_arr_cg[qjd] = UNMETHYL;
                                                    rle_arr_gc[qjd] = IGNORED;
                                                    rle_set = 1;
                                                } else if (qb == 'G') {
                                                    push_char_v(hcg_c, 'C');
                                                    rle_arr_cg[qjd] = METHYLAT;
                                                    rle_arr_gc[qjd] = IGNORED;
                                                    rle_set = 1;
                                                } else {
                                                    push_char_v(hcg_c, 'N');
                                                }
                                            } else if (rb0 != 'C' && rb1 == 'C') { // GCH context
                                                push_int_v(gch_p, (int) rpos+j);
                                                if (qb == 'A') {
                                                    push_char_v(gch_c, 'T');
                                                    rle_arr_cg[qjd] = IGNORED;
                                                    rle_arr_gc[qjd] = SHUT_ACC;
                                                    rle_set = 1;
                                                } else if (qb == 'G') {
                                                    push_char_v(gch_c, 'C');
                                                    rle_arr_cg[qjd] = IGNORED;
                                                    rle_arr_gc[qjd] = OPEN_ACC;
                                                    rle_set = 1;
                                                } else {
                                                    push_char_v(gch_c, 'N');
                                                }
                                            }
                                        }
                                    } else { // bs-seq
                                        rle_arr_gc[qjd] = IGNORED;
                                        if (rb0 == 'C') { // CpG context
                                            // Note: measure G in CpG context
                                            push_int_v(cg_p, (int) rpos+j-1);
                                            if (qb == 'A') {
                                                push_char_v(cg_c, 'T');
                                                rle_arr_cg[qjd] = UNMETHYL;
                                                rle_set = 1;
                                            } else if (qb == 'G') {
                                                push_char_v(cg_c, 'C');
                                                rle_arr_cg[qjd] = METHYLAT;
                                                rle_set = 1;
                                            } else {
                                                push_char_v(cg_c, 'N');
                                            }
                                        }
                                    }
                                }

                                // reference is a C
                                if (!bsstrand && rb == 'C' && rpos+j+1 <= rs->end) {
                                    char rb1 = refcache_getbase_upcase(rs, rpos+j+1); // next base
                                    if (conf->comm.is_nome) { // nome-seq
                                        if (rpos+j-1 >= rs->beg) { // to prevent underflow
                                            char rb0 = refcache_getbase_upcase(rs, rpos+j-1); // previous base
                                            if (rb0 != 'G' && rb1 == 'G') { // HCG context
                                                // measure C in CpG context
                                                push_int_v(hcg_p, (int) rpos+j);
                                                if (qb == 'T') {
                                                    push_char_v(hcg_c, 'T');
                                                    rle_arr_cg[qjd] = UNMETHYL;
                                                    rle_arr_gc[qjd] = IGNORED;
                                                    rle_set = 1;
                                                } else if (qb == 'C') {
                                                    push_char_v(hcg_c, 'C');
                                                    rle_arr_cg[qjd] = METHYLAT;
                                                    rle_arr_gc[qjd] = IGNORED;
                                                    rle_set = 1;
                                                } else {
                                                    push_char_v(hcg_c, 'N');
                                                }
                                            } else if (rb0 == 'G' && rb1 != 'G') { // GCH context
                                                push_int_v(gch_p, (int) rpos+j);
                                                if (qb == 'T') {
                                                    push_char_v(gch_c, 'T');
                                                    rle_arr_cg[qjd] = IGNORED;
                                                    rle_arr_gc[qjd] = SHUT_ACC;
                                                    rle_set = 1;
                                                } else if (qb == 'C') {
                                                    push_char_v(gch_c, 'C');
                                                    rle_arr_cg[qjd] = IGNORED;
                                                    rle_arr_gc[qjd] = OPEN_ACC;
                                                    rle_set = 1;
                                                } else {
                                                    push_char_v(gch_c, 'N');
                                                }
                                            }
                                        }
                                    } else { // bs-seq
                                        rle_arr_gc[qjd] = IGNORED;
                                        if (rb1 == 'G') { // CpG context
                                            push_int_v(cg_p, (int) rpos+j);
                                            if (qb == 'T') {
                                                push_char_v(cg_c, 'T');
                                                rle_arr_cg[qjd] = UNMETHYL;
                                                rle_set = 1;
                                            } else if (qb == 'C') {
                                                push_char_v(cg_c, 'C');
                                                rle_arr_cg[qjd] = METHYLAT;
                                                rle_set = 1;
                                            } else {
                                                push_char_v(cg_c, 'N');
                                            }
                                        }
                                    }
                                }
                            }

                            // Check for SNP
                            uint32_t snp_ind = rpos + j - snp_beg;
                            if (snps && episnp_test(snps, snp_ind)) {
                                // Old formats
                                push_char_v(snp_c, qb);
                                push_int_v(snp_p, rpos+j);

                                // epiBED
                                if (!rle_set) {
                                    rle_arr_cg[qjd] = IGNORED;
                                    rle_arr_gc[qjd] = IGNORED;
                                }

                                if (rle_set && !(episnp_test(meth, snp_ind))) {
                                    rle_arr_cg[qjd] = IGNORED;
                                    rle_arr_gc[qjd] = IGNORED;
                                }
                                if (bsstrand && qb == 'A') {
                                    rle_arr_vr[qjd] = AMBIG_GA;
                                } else if (!bsstrand && qb == 'T') {
                                    rle_arr_vr[qjd] = AMBIG_CT;
                                } else {
                                    rle_arr_vr[qjd] = qb;
                                }

                                rle_set = 1;
                            } else {
                                rle_arr_vr[qjd] = IGNORED;
                                if (!rle_set) {
                                    rle_arr_cg[qjd] = IGNORED;
                                    rle_arr_gc[qjd] = IGNORED;
                                }
                            }

                            // Fill any locations that weren't already filled
                            if (!rle_set) {
                                rle_arr_cg[qjd] = IGNORED;
                                rle_arr_gc[qjd] = IGNORED;
                            } else {
                                rle_set = 0;
                            }
                        }
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CINS:
                        for (j=0; j<oplen; ++j) {
                            qj = qpos + j;
                            qjd = qj + n_deletions;
                            qb = bscall(b, qj);
                            rle_arr_vr[qjd] = tolower(qb);
                            rle_arr_cg[qjd] = SKIP_INS;
                            rle_arr_gc[qjd] = SKIP_INS;

                            if (conf->use_modbam) {
                                // Retrieve next base modification
                                n_mods = bam_mods_at_next_pos(b, mod_state, mod, 1);
                                if (n_mods < 0) {
                                    wzfatal("ERROR: problem encountered retrieving next base modification\n");
                                }
                            }
                        }
                        n_insertions += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CDEL:
                        for (j=0; j<oplen; ++j) {
                            qj = qpos + j;
                            qjd = qj + n_deletions;
                            rle_arr_cg[qjd] = SKIP_DEL;
                            rle_arr_gc[qjd] = SKIP_DEL;
                            rle_arr_vr[qjd] = DELETION;
                        }
                        n_deletions += oplen;
                        rpos += oplen;
                        break;
                    case BAM_CSOFT_CLIP:
                        for (j=0; j<oplen; ++j) {
                            qj = qpos + j;
                            qjd = qj + n_deletions;
                            if (qj <= softclip_start) // We only want to count the softclips when they occur at the start of the read
                                softclip_start++;
                            rle_arr_cg[qjd] = SOFTCLIP;
                            rle_arr_gc[qjd] = SOFTCLIP;
                            rle_arr_vr[qjd] = SOFTCLIP;


                            if (conf->use_modbam) {
                                // Retrieve next base modification
                                n_mods = bam_mods_at_next_pos(b, mod_state, mod, 1);
                                if (n_mods < 0) {
                                    wzfatal("ERROR: problem encountered retrieving next base modification\n");
                                }
                            }
                        }
                        qpos += oplen;
                        break;
                    default:
                        fprintf(stderr, "Unknown cigar %u\n", op);
                        abort();
                }
            }

            // BAM position for reads with softclipping at the beginning is relative the the first Match location
            // in the CIGAR string, so we need to shift backwards the number of softclip positions to get the corresponding
            // start position for the whole read
            uint32_t start = c->pos+1 - softclip_start;
            uint32_t end   = start + c->l_qseq + n_deletions - n_insertions - 1; // -1 adjusts end position due to 1-based counting

            // produce epiread output
            if (conf->epiread_pair) {
                format_epiread_pairwise(
                        &rec.s, chrm, &w, conf,
                        snp_p, hcg_p, gch_p, cg_p,
                        snp_c, hcg_c, gch_c, cg_c);
            }
            if (conf->epiread_old) {
                format_epiread_old(
                        &rec.s, b, bsstrand, chrm, &w, snps, conf,
                        snp_p, hcg_p, gch_p, cg_p,
                        snp_c, hcg_c, gch_c, cg_c);
            }
            if (!conf->epiread_pair && !conf->epiread_old) {
                format_epi_bed(&rec.s, b, bsstrand, chrm, &w, conf, rle_arr_cg, rle_arr_gc, rle_arr_vr, start, end);
            }

            // clean up
            free(rle_arr_gc);
            free(rle_arr_vr);
            free(rle_arr_cg);
            free_int_v(snp_p); free_char_v(snp_c);
            if (conf->comm.is_nome) {
                free_int_v(hcg_p); free_char_v(hcg_c);
                free_int_v(gch_p); free_char_v(gch_c);
            } else {
                free_int_v(cg_p); free_char_v(cg_c);
            }
        }

        // run through cytosines
        rec.block_id = w.block_id;
        // put output string to output queue
        wqueue_put2(record, res->rq, rec);

        if (mod_state) { hts_base_mod_state_free(mod_state); }
        bam_destroy1(b);
        hts_itr_destroy(iter);
        free(meth);
        free(snps);
    }

    free_refcache(rs);
    hts_close(in);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);

    return 0;
}

episnp_chrom1_v *bed_init_episnp(char *snp_bed_fn) {
    episnp_chrom1_v *episnp = init_episnp_chrom1_v(2);
    kstring_t line;
    line.l = line.m = 0; line.s = 0;

    // Read snp bed file
    episnp_chrom1_t *episnp1 = 0;
    char *tok;
    char *ref;
    char *alt;
    double vaf;
    gzFile fh = gzopen(snp_bed_fn, "r");

    if (fh == NULL) {
        free(line.s);
        wzfatal("Could not find SNP BED: %s\n", snp_bed_fn);
    }

    uint8_t first_char = 1;
    while (1) {
        int c = gzgetc(fh);
        if (c < 0 && first_char) {
            free(line.s);
            wzfatal("SNP BED (%s) is empty\n", snp_bed_fn);
        }
        first_char = 0;
        if (c=='\n' || c==EOF) {
            if (strcount_char(line.s, '\t')==8) {
                // Get chromosome
                tok = strtok(line.s, "\t"); // chromosome
                if (!episnp1 || strcmp(episnp1->chrm, tok) != 0)
                    episnp1 = get_n_insert_episnp1(episnp, tok);

                // Adjust number of elements
                episnp1->locs = realloc(episnp1->locs, (episnp1->n+1)*sizeof(uint32_t));
                if (!episnp1->locs) {
                    wzfatal("Failed to properly allocate space for SNP locations\n");
                }
                episnp1->meth = realloc(episnp1->meth, (episnp1->n+1)*sizeof(uint8_t));
                if (!episnp1->meth) {
                    wzfatal("Failed to properly allocate space for SNP methylation calls\n");
                }

                // Get start
                tok = strtok(NULL, "\t"); // start
                ensure_number(tok);
                episnp1->locs[episnp1->n] = atoi(tok)+1;

                // Get ref and alt
                tok = strtok(NULL, "\t"); // end
                tok = strtok(NULL, "\t"); // ref
                ref = tok;
                tok = strtok(NULL, "\t"); // alt
                alt = tok;

                // Get allele frequency
                tok = strtok(NULL, "\t"); // genotype (in biscuit)
                tok = strtok(NULL, "\t"); // allele support (in biscuit)
                tok = strtok(NULL, "\t"); // allele depth for allele frequency
                tok = strtok(NULL, "\t"); // allele frequency
                ensure_number(tok);
                vaf = atof(tok);

                uint8_t meth_callable = 0;
                if (strcmp(ref, "C") == 0) {
                    if ((strcmp(alt, "T") != 0) || ((strcmp(alt, "T") == 0) && (vaf < 0.05))) {
                        meth_callable = 1;
                    }
                }
                if (strcmp(ref, "G") == 0) {
                    if ((strcmp(alt, "A") != 0) || ((strcmp(alt, "A") == 0) && (vaf < 0.05))) {
                        meth_callable = 1;
                    }
                }

                episnp1->meth[episnp1->n] = meth_callable;
                episnp1->n++;
            }

            line.l = 0;
            if (c==EOF) {
                break;
            }
        } else {
            kputc(c, &line);
        }
    }
    free(line.s);

    gzclose(fh);

    return episnp;
}

void epiread_conf_init(epiread_conf_t *conf) {
    conf->comm = bisc_common_init();
    conf->bt = bisc_threads_init();
    conf->filt = meth_filter_init();

    conf->epiread_reg_start = 0;
    conf->epiread_reg_end = 0;
    conf->modbam_prob = 0.9;
    conf->filter_empty_epiread = 1;
    conf->is_long_read = 0;
    conf->epiread_old = 0;
    conf->epiread_pair = 0;
    conf->print_all_locations = 0;
    conf->use_modbam = 0;
}

static int usage() {
    epiread_conf_t conf;
    epiread_conf_init(&conf);

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit epiread [options] <ref.fa> <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -B STR    Bed input for SNP display in epiread output\n");
    fprintf(stderr, "    -g STR    Region (optional, will process the whole bam if not specified)\n");
    fprintf(stderr, "    -s STR    Step of window dispatching [%d]\n", conf.bt.step);
    fprintf(stderr, "    -@ INT    Number of threads [%d]\n", conf.bt.n_threads);
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o STR    Output file [stdout]\n");
    fprintf(stderr, "    -N        NOMe-seq mode [off]\n");
    fprintf(stderr, "    -L        Data is from long read sequencing [off]\n");
    fprintf(stderr, "    -M        BAM file has modBAM tags (MM/ML) [off]\n");
    fprintf(stderr, "    -P        Pairwise mode [off]\n");
    fprintf(stderr, "    -O        Old BISCUIT epiread format, not compatible with -P [off]\n");
    fprintf(stderr, "    -A        Print all CpG and SNP locations in location column, ignored if -O not given [off]\n");
    fprintf(stderr, "    -v        Verbose (print additional info for diagnostics) [off]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -b INT    Minimum base quality [%u]\n", conf.filt.min_base_qual);
    fprintf(stderr, "    -m INT    Minimum mapping quality [%u]\n", conf.filt.min_mapq);
    fprintf(stderr, "    -a INT    Minimum alignment score (from AS-tag) [%u]\n", conf.filt.min_score);
    fprintf(stderr, "    -t INT    Max cytosine retention in a read [%u]\n", conf.filt.max_retention);
    fprintf(stderr, "    -l INT    Minimum read length [%u]\n", conf.filt.min_read_len);
    fprintf(stderr, "    -5 INT    Minimum distance to 5' end of a read [%u]\n", conf.filt.min_dist_end_5p);
    fprintf(stderr, "    -3 INT    Minimum distance to 3' end of a read [%u]\n", conf.filt.min_dist_end_3p);
    fprintf(stderr, "    -E        NO filtering of empty epireads\n");
    //fprintf(stderr, "    -c        NO filtering secondary mapping\n");
    fprintf(stderr, "    -d        Double count cytosines in overlapping mate reads (avoided\n");
    fprintf(stderr, "                  by default)\n");
    fprintf(stderr, "    -u        NO filtering of duplicate\n");
    fprintf(stderr, "    -p        NO filtering of improper pair\n");
    fprintf(stderr, "    -n INT    Maximum NM tag [%d]\n", conf.filt.max_nm);
    fprintf(stderr, "    -y FLT    Minimum probability a modification is correct (0.0 - 1.0) [%f]\n", conf.modbam_prob);
    fprintf(stderr, "    -h        This help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note, the -O (old epiread format) and -P (pairwise format for biscuit asm) are not guaranteed\n");
    fprintf(stderr, "    to match output from biscuit pileup. These file formats have been left in for legacy purposes.\n");
    fprintf(stderr, "    Default output with an unfiltered BISCUIT SNP BED file (biscuit pileup ...\n");
    fprintf(stderr, "    -> biscuit vcf2bed -t snp ...) should have the same results in the epiBED as in pileup.\n");
    fprintf(stderr, "\n");

    return 1;
}

int main_epiread(int argc, char *argv[]) {
    int c;
    char *reg = 0;
    char *outfn = 0;
    char *statsfn = 0;
    char *snp_bed_fn = 0;
    epiread_conf_t conf;
    epiread_conf_init(&conf);

    if (argc<2) return usage();
    while ((c=getopt(argc, argv, ":@:B:o:g:s:t:l:5:3:n:b:m:a:y:AMNLEcduOPpvh"))>=0) {
        switch (c) {
            case 'B': snp_bed_fn = optarg; break;
            case 'o': outfn = optarg; break;
            case 'g': reg = optarg; break;
            case '@': conf.bt.n_threads = atoi(optarg); break;
            case 's': conf.bt.step = atoi(optarg); break;
            case 't': conf.filt.max_retention = atoi(optarg); break;
            case 'l': conf.filt.min_read_len = atoi(optarg); break;
            case '5': conf.filt.min_dist_end_5p = atoi(optarg); break;
            case '3': conf.filt.min_dist_end_3p = atoi(optarg); break;
            case 'n': conf.filt.max_nm = atoi(optarg); break;
            case 'b': conf.filt.min_base_qual = atoi(optarg); break;
            case 'm': conf.filt.min_mapq = atoi(optarg); break;
            case 'a': conf.filt.min_score = atoi(optarg); break;
            case 'O': conf.epiread_old = 1; break;
            case 'A': conf.print_all_locations = 1; break;
            case 'N': conf.comm.is_nome = 1; break;
            case 'L': conf.is_long_read = 1; break;
            case 'M': conf.use_modbam = 1; break;
            case 'y': conf.modbam_prob = atof(optarg); break;
            case 'E': conf.filter_empty_epiread = 0; break;
            case 'c': conf.filt.filter_secondary = 0; break;
            case 'd': conf.filt.filter_doublecnt = 0; break;
            case 'u': conf.filt.filter_duplicate = 0; break;
            case 'p': conf.filt.filter_ppair = 0; break;
            case 'P': conf.epiread_pair = 1; break;
            case 'v': conf.comm.verbose = 1; break;
            case 'h': return usage();
            case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt); break;
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt); break;
            default: return usage();
        }
    }

    if (conf.epiread_old && conf.epiread_pair) {
        usage();
        wzfatal("Cannot run with both pairwise and old epiread format set.\n");
    }

    if (conf.modbam_prob < 0.0 || conf.modbam_prob > 1.0) {
        usage();
        wzfatal("Minimum modification probability must be between 0.0 and 1.0\n");
    }

    if (optind + 2 > argc) {
        usage();
        wzfatal("Reference or bam input is missing\n");
    }
    char *reffn = argv[optind++];
    char *infn = argv[optind++];

    episnp_chrom1_v *episnp = snp_bed_fn ?
        bed_init_episnp(snp_bed_fn) : NULL;

    wqueue_t(window) *wq = wqueue_init(window, 100000);
    pthread_t *processors = calloc(conf.bt.n_threads, sizeof(pthread_t));
    result_t *results = calloc(conf.bt.n_threads, sizeof(result_t));
    int i; unsigned j;
    htsFile *in = hts_open(infn, "rb");
    if (in == NULL) {
        wzfatal("%s unable to be opened\n", infn);
    }
    bam_hdr_t *header = sam_hdr_read(in);

    // sort sequence name by alphabetic order, chr1, chr10, chr11 ...
    target_v *targets = init_target_v(50);
    target_t *t;
    for (i=0; i<header->n_targets; ++i) {
        t = next_ref_target_v(targets);
        t->tid = i;
        t->name = header->target_name[i];
        t->len = header->target_len[i];
    }

    qsort(targets->buffer, targets->size,
            sizeof(target_t), compare_targets);

    // setup writer
    pthread_t writer;
    writer_conf_t writer_conf = {
        .q = wqueue_init(record, 100000),
        .outfn = outfn,
        .statsfn = statsfn,
        .header = 0,
        .targets = targets,
        .conf = &conf,
    };
    pthread_create(&writer, NULL, epiread_write_func, &writer_conf);
    for (i=0; i<conf.bt.n_threads; ++i) {
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        results[i].snp = episnp;
        results[i].ref_fn = reffn;
        results[i].bam_fn = infn;
        results[i].conf = &conf;
        pthread_create(&processors[i], NULL, process_func, &results[i]);
    }

    window_t w; memset(&w, 0, sizeof(window_t));
    uint32_t wbeg;
    int64_t block_id=0;

    // process bam
    if (reg) { // regional
        int tid;
        uint32_t beg, end;
        biscuit_parse_region(reg, header, &tid, (int*) &beg, (int*) &end);
        // chromosome are assumed to be less than 2**29
        beg++; // shift beg from 0-based to 1-based
        if (beg<=0) beg = 1;
        if (end>header->target_len[tid]) end = header->target_len[tid];

        // Save the start and end location when a region is called to catch data that falls outside these regions
        // This is correct behavior based on how samtools finds reads that overlap a specified region
        conf.epiread_reg_start = beg;
        conf.epiread_reg_end   = end;

        for (wbeg = beg; wbeg < end; wbeg += conf.bt.step, block_id++) {
            w.tid = tid;
            w.block_id = block_id;
            w.beg = wbeg;
            w.end = wbeg + conf.bt.step;
            if (w.end > end) w.end = end;
            wqueue_put(window, wq, &w);
        }
    } else { // entire bam
        for (j=0; j<targets->size; ++j) {
            t = ref_target_v(targets, j);
            for (wbeg = 1; wbeg < t->len; wbeg += conf.bt.step, block_id++) {
                w.tid = t->tid;
                w.block_id = block_id;
                w.beg = wbeg;
                w.end = wbeg+conf.bt.step;
                if (w.end > t->len) w.end = t->len;
                wqueue_put(window, wq, &w);
            }
        }
    }
    for (i=0; i<conf.bt.n_threads; ++i) {
        w.tid = -1;
        wqueue_put(window, wq, &w);
    }

    for (i=0; i<conf.bt.n_threads; ++i) {
        pthread_join(processors[i], NULL);
    }

    record_t rec = { .block_id = RECORD_QUEUE_END };
    wqueue_put2(record, writer_conf.q, rec);
    pthread_join(writer, NULL);
    wqueue_destroy(record, writer_conf.q);

    free_target_v(targets);
    free(results);
    free(processors);
    wqueue_destroy(window, wq);
    hts_close(in);
    bam_hdr_destroy(header);

    if (episnp)
        destroy_episnp(episnp);

    return 0;
}
