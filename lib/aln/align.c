/* entry point for alignment
 *
 * Newly added copyright in 2022
 * Copyright (c) 2022-2023 Jacob.Morrison@vai.org
 *
 * The MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
#include "wzmisc.h"
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
  kseq_t *ks, *ks2;
  mem_opt_t *opt;
  mem_pestat_t *pes0;
  int64_t n_processed;
  int copy_comment, actual_chunk_size;
  // for debug
  char *_seq1, *_seq2;
  int __processed;
  bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
  ktp_aux_t *aux;
  int n_seqs;
  bseq1_t *seqs;
} ktp_data_t;

static void *process(void *shared, int step, void *_data) {
  ktp_aux_t *aux = (ktp_aux_t*)shared;
  ktp_data_t *data = (ktp_data_t*)_data;
  int i;
  if (step == 0) {
    int64_t size = 0;
    ktp_data_t *ret = calloc(1, sizeof(ktp_data_t));
    if (aux->_seq1) {  // prompt supplied input, for debug
      if (aux->__processed) return 0;
      ret->seqs = bis_create_bseq1(aux->_seq1, aux->_seq2, &ret->n_seqs);
      if (aux->_seq2)  aux->opt->flag |= MEM_F_PE;
      aux->__processed = 1;
    } else { // read from file
      ret->seqs = bis_bseq_read(aux->actual_chunk_size, aux->opt->has_bc, &ret->n_seqs, aux->ks, aux->ks2);
      if (ret->seqs == 0) {
        free(ret);
        return 0;
      }

      if (!aux->copy_comment)
        for (i = 0; i < ret->n_seqs; ++i) {
          free(ret->seqs[i].comment);
          ret->seqs[i].comment = 0;
        }
    }

    for (i = 0; i < ret->n_seqs; ++i) {
      ret->seqs[i].sam = 0;
      size += ret->seqs[i].l_seq;
    }
    
    if (bwa_verbose >= 3)
      fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
    
    return ret;
  } else if (step == 1) {
    const mem_opt_t *opt = aux->opt;
    const bwaidx_t *idx = aux->idx;

    /* interleaved input */
    if (opt->flag & MEM_F_SMARTPE) {
      bseq1_t *sep[2];
      int n_sep[2];
      mem_opt_t tmp_opt = *opt;

      bseq_classify(data->n_seqs, data->seqs, n_sep, sep);

      if (n_sep[0]) {           // single-end
        tmp_opt.flag &= ~MEM_F_PE;
        mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, n_sep[0], sep[0], 0);
        for (i = 0; i < n_sep[0]; ++i)
          data->seqs[sep[0][i].id].sam = sep[0][i].sam;
      }

      if (n_sep[1]) {           // paired-end
        tmp_opt.flag |= MEM_F_PE;
        mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0);
        for (i = 0; i < n_sep[1]; ++i)
          data->seqs[sep[1][i].id].sam = sep[1][i].sam;
      }

      // clean up
      if (n_sep[0]) {
          for (i = 0; i < n_sep[0]; ++i) {
              /* bisulfite free, the pointers can be NULL */
              free(sep[0][i].bisseq[0]);
              free(sep[0][i].bisseq[1]);
          }
      }
      if (n_sep[1]) {
          for (i = 0; i < n_sep[1]; ++i) {
              /* bisulfite free, the pointers can be NULL */
              free(sep[1][i].bisseq[0]);
              free(sep[1][i].bisseq[1]);
          }
      }
      free(sep[0]); free(sep[1]);
    } else {
      mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seqs, aux->pes0);
    }

    aux->n_processed += data->n_seqs;

    return data;
  } else if (step == 2) {
    for (i = 0; i < data->n_seqs; ++i) {
      if (data->seqs[i].sam) err_fputs(data->seqs[i].sam, stdout);
      free(data->seqs[i].name); free(data->seqs[i].comment);
      free(data->seqs[i].seq0); free(data->seqs[i].qual); free(data->seqs[i].sam);
      /* bisulfite free, the pointers can be NULL */
      free(data->seqs[i].bisseq[0]);
      free(data->seqs[i].bisseq[1]);
    }
    free(data->seqs); free(data);
    return 0;
  }
  return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0) {
  if (opt0->a) { // matching score is changed
    if (!opt0->b) opt->b *= opt->a;
    if (!opt0->T) opt->T *= opt->a;
    if (!opt0->o_del) opt->o_del *= opt->a;
    if (!opt0->e_del) opt->e_del *= opt->a;
    if (!opt0->o_ins) opt->o_ins *= opt->a;
    if (!opt0->e_ins) opt->e_ins *= opt->a;
    if (!opt0->zdrop) opt->zdrop *= opt->a;
    if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
    if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
    if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
  }
}

static void infer_alt_chromosomes(bntseq_t *bns) {
  // logic: if the chr1,...chr22,chrX,chrY,chrM,
  // then mark everything with chrUn and _random and _hap as alt
  int i;
  for (i=0; i<bns->n_seqs; ++i) if (bns->anns[i].is_alt) break;
  if (i<bns->n_seqs) return;    // skip if alt info is present
  
  int found[25]; memset(found, 0, 25*sizeof(int));
  for (i=0; i<bns->n_seqs; ++i) {
    if (strncmp(bns->anns[i].name, "chr", 3) == 0) {
      if (strlen(bns->anns[i].name) == 4) {
        if (toupper(bns->anns[i].name[3]) == 'X') found[22] = 1;
        else if (toupper(bns->anns[i].name[3]) == 'Y') found[23] = 1;
        else if (toupper(bns->anns[i].name[3]) == 'M') found[24] = 1;
        else if (isdigit(bns->anns[i].name[3])) {
          int n = bns->anns[i].name[3]-'0';
          if (n>0 && n<=22) found[n-1] = 1;
        }
      } else if (strlen(bns->anns[i].name) == 5 &&
                 isdigit(bns->anns[i].name[3]) && isdigit(bns->anns[i].name[4])) {
        int n = atoi(bns->anns[i].name+3);
        if (n>0 && n<=22) found[n-1] = 1;
      }
    }
  }

  int n;
  for (i=n=0; i<25; ++i) if (found[i]) ++n;
  if (n < 20) return; // should take care of human and mouse
  
  for (i=0; i<bns->n_seqs; ++i) {
    if (strncmp(bns->anns[i].name, "chrUn", 5)==0 || 
        strstr(bns->anns[i].name, "_random") || 
        strstr(bns->anns[i].name, "_hap") ||
        strstr(bns->anns[i].name, "_alt")) {
      bns->anns[i].is_alt = 1;
      if (bwa_verbose >= 4)
        fprintf(stderr, "[M:%s] Set %s as ALT.\n", __func__, bns->anns[i].name);
    }
  }
}

int usage(mem_opt_t *opt) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit align [options] <fai-index base> <in1.fq> [in2.fq]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Algorithm options:\n");
    fprintf(stderr, "    -@ INT          Number of threads [%d]\n", opt->n_threads);
    fprintf(stderr, "    -b INT          For PE, map read 1 to parent, read 2 to daughter\n");
    fprintf(stderr, "                        (1: directional library) or read1 and read 2 to both\n");
    fprintf(stderr, "                        (0: non-directional library [default]). For SE, parent\n");
    fprintf(stderr, "                        (1), daughter (3), or both (0 [default]). Note, parent\n");
    fprintf(stderr, "                        is the bisulfite treated strand and the daughter is\n");
    fprintf(stderr, "                        the complement strand.\n");
    fprintf(stderr, "    -f INT          1: BSW strand; 3: BSC strand; 0: both. Note, libraries\n");
    fprintf(stderr, "                        targeting either BSW or BSC are unseen so far! [0]\n");
    fprintf(stderr, "    -k INT          Minimum seed length [%d]\n", opt->min_seed_len);
    fprintf(stderr, "    -w INT          Band width for banded alignment [%d]\n", opt->w);
    fprintf(stderr, "    -d INT          Off-diagonal X-dropoff [%d]\n", opt->zdrop);
    fprintf(stderr, "    -r FLOAT        Look for internal seeds inside a seed longer than\n");
    fprintf(stderr, "                        {-k}*FLOAT [%g]\n", opt->split_factor);
    fprintf(stderr, "    -y INT          Seed occurrence for the 3rd round of seeding [%ld]\n", (long)opt->max_mem_intv);
    /* fprintf(stderr, "    -s INT          Look for internal seeds inside a seed with less than\n"); *
     * fprintf(stderr, "                        INT occ [%d]\n", opt->split_width); */
    fprintf(stderr, "    -J STR          Adaptor of read 1 (fastq direction)\n");
    fprintf(stderr, "    -K STR          Adaptor of read 2 (fastq direction)\n");
    fprintf(stderr, "    -z INT          Minimum base quality to keep from both ends of reads [%d]\n", opt->min_base_qual);
    fprintf(stderr, "    -5 INT          Number of extra bases to clip from 5'-end [%d]\n", opt->clip5);
    fprintf(stderr, "    -3 INT          Number of extra bases to clip from 3'-end [%d]\n", opt->clip3);
    fprintf(stderr, "    -c INT          Skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
    fprintf(stderr, "    -D FLOAT        Drop chains shorter than FLOAT fraction of the longest\n");
    fprintf(stderr, "                        overlapping chain [%.2f]\n", opt->drop_ratio);
    fprintf(stderr, "    -W INT          Discard a chain if seeded bases shorter than INT [0]\n");
    fprintf(stderr, "    -m INT          Perform at most INT rounds of mate rescues for a read [%d]\n", opt->max_matesw);
    fprintf(stderr, "    -S              Skip mate rescue\n");
    fprintf(stderr, "    -P              Skip pairing - mate rescue performed unless -S also given\n");
    fprintf(stderr, "    -e              Discard full-length exact matches\n");
    fprintf(stderr, "    -9              Extract barcode from read comment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Scoring options:\n");
    fprintf(stderr, "    -A INT          Score for a sequence match, scales options -TdBOELU unless\n");
    fprintf(stderr, "                        overridden [%d]\n", opt->a);
    fprintf(stderr, "    -B INT          Penalty for a mismatch [%d]\n", opt->b);
    fprintf(stderr, "    -O INT[,INT]    Gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
    fprintf(stderr, "    -E INT[,INT]    Gap extension penalty; a gap of size, g, has penalty\n");
    fprintf(stderr, "                        {-O} + {-E}*g [%d,%d]\n", opt->e_del, opt->e_ins);
    fprintf(stderr, "    -L INT[,INT]    Penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
    fprintf(stderr, "    -U INT          Penalty for an unpaired read pair [%d]\n", opt->pen_unpaired);
    /* fprintf(stderr, "    -x STR          Read type. Setting -x changes multiple parameters unless\n"); *
     * fprintf(stderr, "                        overriden\n"); *
     * fprintf(stderr, "                        pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0\n"); *
     * fprintf(stderr, "                                (PacBio reads to ref)\n"); *
     * fprintf(stderr, "                        ont2d:  -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0\n"); *
     * fprintf(stderr, "                                (Oxford Nanopore 2D-reads to ref)\n"); *
     * fprintf(stderr, "                        intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n"); *
     * fprintf(stderr, "                        pbread: -k13 -W40 -c1000 -r10 -A1 -B1 -O1 -E1 -N25\n"); *
     * fprintf(stderr, "                                -FeaD.001\n"); */
    fprintf(stderr, "\n");
    fprintf(stderr, "Input/output options:\n");
    fprintf(stderr, "    -1 STR          Align a single read STR\n");
    fprintf(stderr, "    -2 STR          Align a read STR paired with -1 read\n");
    fprintf(stderr, "    -i              Turn off autoinference of ALT chromosomes\n");
    fprintf(stderr, "    -p              Smart pairing (ignores in2.fq)\n");
    fprintf(stderr, "    -R STR          Read group header line (such as '@RG\\tID:foo\\tSM:bar')\n");
    fprintf(stderr, "    -F              Suppress SAM header output\n");
    fprintf(stderr, "    -H STR/FILE     Insert STR to header if it starts with @ or insert lines\n");
    fprintf(stderr, "                        in FILE\n");
    fprintf(stderr, "    -j              Treat ALT contigs as part of the primary assembly (i.e.\n");
    fprintf(stderr, "                        ignore <fai-index base>.alt file)\n");
    fprintf(stderr, "    -q              Do not modify mapQ of supplementary alignments\n");
    fprintf(stderr, "    -T INT          Minimum score to output [%d]\n", opt->T);
    /* if there are <INT> hits with score >80% of the max score, output all in XA *
     *     max_XA_hits     - maximum number of hits on primary chromosomes *
     *     max_XA_hits_alt - maximum number of hits on alt chromosomes */
    fprintf(stderr, "    -g INT[,INT]    Maximum number of hits output in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
    fprintf(stderr, "    -a              Output all alignments for SE or unpaired PE\n");
    fprintf(stderr, "    -C              Append FASTA/FASTQ comment to SAM output\n");
    fprintf(stderr, "    -V              Output the reference FASTA header in the XR tag\n");
    fprintf(stderr, "    -Y              Use soft clipping for supplementary alignments\n");
    fprintf(stderr, "    -M              Mark shorter split hits as secondary\n");
    fprintf(stderr, "    -I FLOAT[,FLOAT[,INT[,INT]]]\n");
    fprintf(stderr, "                    Specify the mean, standard deviation (10%% of the mean\n");
    fprintf(stderr, "                        if absent), maximum (4 sigma from the mean if absent)\n");
    fprintf(stderr, "                        and minimum of insert size distribution. FR orientation\n");
    fprintf(stderr, "                        only [inferred]\n");
    fprintf(stderr, "    -v INT          Verbosity level: \n");
    fprintf(stderr, "                        1: error, 2: warning, 3: message, 4+: debugging [%d]\n", bwa_verbose);
    fprintf(stderr, "    -h              This help\n");
    fprintf(stderr, "\n");

    free(opt);
    return 1;
}

/* the old main_mem */
int main_align(int argc, char *argv[]) {
  mem_opt_t *opt, opt0;
  int fd, fd2, i, c, ignore_alt = 0, no_mt_io = 0;
  int fixed_chunk_size = -1;
  char *p, *rg_line = 0, *hdr_line = 0;
  const char *mode = 0;
  //mem_pestat_t pes[4];
  ktp_aux_t aux;

  memset(&aux, 0, sizeof(ktp_aux_t));
  //memset(pes, 0, 4 * sizeof(mem_pestat_t));
  //for (i = 0; i < 4; ++i) pes[i].failed = 1;


  aux.pes0 = NULL;
  aux.opt = opt = mem_opt_init();
  opt->flag |= MEM_F_NO_MULTI;  /* WZBS */
  memset(&opt0, 0, sizeof(mem_opt_t));
  int auto_infer_alt_chrom = 1;
  if (argc < 2) return usage(opt);
  while ((c = getopt(argc, argv, ":@:1:2:3:5:9ab:c:d:ef:g:hijk:m:pqr:s:v:w:x:y:z:A:B:CD:E:FG:H:I:J:K:L:MN:O:PQ:R:ST:U:VW:X:Y")) >= 0) {
      if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
      else if (c == '1') aux._seq1 = strdup(optarg);
      else if (c == '2') aux._seq2 = strdup(optarg);
      else if (c == 'x') mode = optarg;
      else if (c == 'b') opt->parent = atoi(optarg);   /* targeting parent or daughter */
      else if (c == 'f') opt->bsstrand = atoi(optarg); /* targeting BSW or BSC */
      else if (c == 'i') auto_infer_alt_chrom = 0; // turn off auto-inference of alt-chromosomes
      else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
      else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
      else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
      else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
      else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
      else if (c == '@') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
      else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
      else if (c == 'a') opt->flag |= MEM_F_ALL;
      else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
      else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
      else if (c == 'M') opt->flag |= MEM_F_NO_MULTI; // mark shorter split as secondary
      else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
      else if (c == 'e') opt->flag |= MEM_F_SELF_OVLP;
      else if (c == 'F') opt->flag |= MEM_F_ALN_REG;
      else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
      else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
      else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
      else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
      else if (c == 'v') bwa_verbose = atoi(optarg);
      else if (c == 'j') ignore_alt = 1;
      else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
      else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
      else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
      else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
      else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
      else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
      else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
      else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
      else if (c == 'C') aux.copy_comment = 1;
      /* else if (c == 'K') fixed_chunk_size = atoi(optarg); -K now reads adapter sequences for read 2*/
      else if (c == 'J') {
          opt->l_adaptor1 = strlen(optarg);
          opt->adaptor1 = calloc(opt->l_adaptor1, sizeof(uint8_t));
          for (i=0; i<opt->l_adaptor1; ++i)
              opt->adaptor1[i] = nst_nt4_table[(int)optarg[i]];
      } else if (c == 'K') {
          opt->l_adaptor2 = strlen(optarg);
          opt->adaptor2 = calloc(opt->l_adaptor2, sizeof(uint8_t));
          for (i=0; i<opt->l_adaptor2; ++i)
              opt->adaptor2[i] = nst_nt4_table[(int)optarg[i]];
      } else if (c == 'z') opt->min_base_qual = atoi(optarg);
      else if (c == '5') opt->clip5 = atoi(optarg);
      else if (c == '3') opt->clip3 = atoi(optarg);
      else if (c == '9') opt->has_bc = 1;
      else if (c == 'X') opt->mask_level = atof(optarg);
      else if (c == 'g') {
          opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
          opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
          if (*p != 0 && ispunct(*p) && isdigit(p[1]))
              opt->max_XA_hits_alt = strtol(p+1, &p, 10);
      }
      else if (c == 'Q') {
          opt0.mapQ_coef_len = 1;
          opt->mapQ_coef_len = atoi(optarg);
          opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
      } else if (c == 'O') {
          opt0.o_del = opt0.o_ins = 1;
          opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
          if (*p != 0 && ispunct(*p) && isdigit(p[1]))
              opt->o_ins = strtol(p+1, &p, 10);
      } else if (c == 'E') {
          opt0.e_del = opt0.e_ins = 1;
          opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
          if (*p != 0 && ispunct(*p) && isdigit(p[1]))
              opt->e_ins = strtol(p+1, &p, 10);
      } else if (c == 'L') {
          opt0.pen_clip5 = opt0.pen_clip3 = 1;
          opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
          if (*p != 0 && ispunct(*p) && isdigit(p[1]))
              opt->pen_clip3 = strtol(p+1, &p, 10);
      } else if (c == 'R') {
          if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
      } else if (c == 'H') { // header to be copied from the given file
          if (optarg[0] != '@') {
              FILE *fp;
              if ((fp = fopen(optarg, "r")) != 0) {
                  char *buf;
                  buf = calloc(1, 0x10000);
                  while (fgets(buf, 0xffff, fp)) {
                      i = strlen(buf);
                      assert(buf[i-1] == '\n'); // a long line
                      buf[i-1] = 0;
                      hdr_line = bwa_insert_header(buf, hdr_line);
                  }
                  free(buf);
                  fclose(fp);
              }
          } else hdr_line = bwa_insert_header(optarg, hdr_line);
      } else if (c == 'I') { // specify the insert size distribution
          mem_pestat_t *pes = calloc(1, sizeof(mem_pestat_t));
          pes->avg = strtod(optarg, &p);
          pes->avg = strtod(optarg, &p);
          pes->std = pes->avg * .1;
          if (*p != 0 && ispunct(*p) && isdigit(p[1]))
              pes->std = strtod(p+1, &p);
          pes->high = (int)(pes->avg + 4. * pes->std + .499);
          pes->low  = (int)(pes->avg - 4. * pes->std + .499);
          //if (pes->low < 1) pes->low = 1;
          if (*p != 0 && ispunct(*p) && isdigit(p[1]))
              pes->high = (int)(strtod(p+1, &p) + .499);
          if (*p != 0 && ispunct(*p) && isdigit(p[1]))
              pes->low  = (int)(strtod(p+1, &p) + .499);
          if (bwa_verbose >= 3)
              fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
                      __func__, pes->avg, pes->std, pes->high, pes->low);
          aux.pes0 = pes;
      } else if (c == 'h') {
          return usage(opt);
      } else if (c == ':') {
          usage(opt); wzfatal("Option needs an argument: -%c\n", optopt);
      } else if (c == '?') {
          usage(opt); wzfatal("Unrecognized option: -%c\n", optopt);
      } else {
          return usage(opt);
      }
  }

  if (rg_line) {
    hdr_line = bwa_insert_header(rg_line, hdr_line);
    free(rg_line);
  }

  if (opt->n_threads < 1) opt->n_threads = 1;
  if ((optind + 1 >= argc || optind + 3 < argc) && !aux._seq1) {
      usage(opt);
      wzfatal("Missing fai-index base or FASTQ file\n");
  }

  // different mode for different library prep
  if (mode) {
    if (strcmp(mode, "intractg") == 0) {
      if (!opt0.o_del) opt->o_del = 16;
      if (!opt0.o_ins) opt->o_ins = 16;
      if (!opt0.b) opt->b = 9;
      if (!opt0.pen_clip5) opt->pen_clip5 = 5;
      if (!opt0.pen_clip3) opt->pen_clip3 = 5;
    } else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "pbread") == 0 || strcmp(mode, "ont2d") == 0) {
      if (!opt0.o_del) opt->o_del = 1;
      if (!opt0.e_del) opt->e_del = 1;
      if (!opt0.o_ins) opt->o_ins = 1;
      if (!opt0.e_ins) opt->e_ins = 1;
      if (!opt0.b) opt->b = 1;
      if (opt0.split_factor == 0.) opt->split_factor = 10.;
      if (strcmp(mode, "pbread") == 0) { // pacbio read-to-read setting; NOT working well!
        opt->flag |= MEM_F_ALL | MEM_F_SELF_OVLP | MEM_F_ALN_REG;
        if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
        if (!opt0.max_occ) opt->max_occ = 1000;
        if (!opt0.min_seed_len) opt->min_seed_len = 13;
        if (!opt0.max_chain_extend) opt->max_chain_extend = 25;
        if (opt0.drop_ratio == 0.) opt->drop_ratio = .001;
      } else if (strcmp(mode, "ont2d") == 0) {
        if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
        if (!opt0.min_seed_len) opt->min_seed_len = 14;
        if (!opt0.pen_clip5) opt->pen_clip5 = 0;
        if (!opt0.pen_clip3) opt->pen_clip3 = 0;
      } else {
        if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
        if (!opt0.min_seed_len) opt->min_seed_len = 17;
        if (!opt0.pen_clip5) opt->pen_clip5 = 0;
        if (!opt0.pen_clip3) opt->pen_clip3 = 0;
      }
    } else {
      fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
      return 1; // FIXME memory leak
    }
  } else update_a(opt, &opt0);

  /* setup scoring */
  bwa_fill_scmat(opt->a, opt->b, opt->mat);
  bwa_fill_scmat_ct(opt->a, opt->b, opt->ctmat);
  bwa_fill_scmat_ga(opt->a, opt->b, opt->gamat);

  /* load bwt index */
  if (optind >= argc) {
      usage(opt);
      wzfatal("Missing fai-index base\n");
  }
  aux.idx = bwa_idx_load_from_shm(argv[optind]);
  if (aux.idx == 0) {
    if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
  } else if (bwa_verbose >= 3)
    fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);

  // infer alternative chromosomes from name
  if (auto_infer_alt_chrom) infer_alt_chromosomes(aux.idx->bns);
  
  if (ignore_alt)
    for (i = 0; i < aux.idx->bns->n_seqs; ++i)
      aux.idx->bns->anns[i].is_alt = 0;

  gzFile fp, fp2 = 0;
  void *ko = 0, *ko2 = 0;
  if  (!aux._seq1) {
    /* setup fastq input */
    ko = kopen(argv[optind + 1], &fd);
    if (ko == 0) {
      if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
      return 1;
    }
    fp = gzdopen(fd, "r");
    aux.ks = kseq_init(fp);
    if (optind + 2 < argc) {
      if (opt->flag&MEM_F_PE) {
        if (bwa_verbose >= 2)
          fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
      } else {
        ko2 = kopen(argv[optind + 2], &fd2);
        if (ko2 == 0) {
          if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
          return 1;
        }
        fp2 = gzdopen(fd2, "r");
        aux.ks2 = kseq_init(fp2);
        opt->flag |= MEM_F_PE;
      }
    }
  } else {fp = fp2 = NULL; ko = ko2 = NULL;}

  /* print header */
  if (!(opt->flag & MEM_F_ALN_REG))
    bwa_print_sam_hdr(aux.idx->bns, hdr_line);
  aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
  kt_pipeline(no_mt_io? 1 : 2, process, &aux, 3);
  free(hdr_line);
  free(opt->adaptor1); free(opt->adaptor2);
  free(opt);
  bwa_idx_destroy(aux.idx);
  kseq_destroy(aux.ks);
  if (fp) {
    err_gzclose(fp);
    kclose(ko);
    if (ko) free(ko);           /* kclose doesn't do that */
  }
  free(aux.pes0);
  if (aux.ks2) {
    kseq_destroy(aux.ks2);
    if (fp2) {
      err_gzclose(fp2);
      kclose(ko2);
      if (ko2) free(ko2);       /* kclose doesn't do that */
    }
  }
  return 0;
}
