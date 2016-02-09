/**
 * convert bam to epiread format with supplied SNP bed file 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Wanding.Zhou@vai.org
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
**/
#include "pileup.h"
#include "wstr.h"

typedef struct episnp_chrom1_t {
  char *chrm;
  size_t n;
  uint32_t *locs;
} episnp_chrom1_t;

DEFINE_VECTOR(episnp_chrom1_v, episnp_chrom1_t)

void destroy_episnp(episnp_chrom1_v *episnp) {
  uint32_t i;
  for (i=0; i<episnp->size; ++i) {
    episnp_chrom1_t *e = ref_episnp_chrom1_v(episnp, i);
    free(e->locs);
  }
  free_episnp_chrom1_v(episnp);
}

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
    episnp1->n = 0;
  }
  return episnp1;
}

typedef struct {
  char *bam_fn;                 /* on stack */
  char *ref_fn;                 /* on stack */
  episnp_chrom1_v *snp;
  wqueue_t(window) *q;
  wqueue_t(record) *rq;
  conf_t *conf;
} result_t;

static void *epiread_write_func(void *data) {
  writer_conf_t *c = (writer_conf_t*) data;

  FILE *out;
  if (c->outfn) out=fopen(c->outfn, "w");
  else out=stdout;

  int64_t next_block = 0;
  record_v *records = init_record_v(20);

  while (1) {
    record_t rec;
    wqueue_get(record, c->q, &rec);
    if(rec.block_id == RECORD_QUEUE_END) break;
    if (rec.block_id == next_block) {
      do {
        if (rec.s.s)
	  fputs(rec.s.s, out);
        free(rec.s.s);

        /* get next block from shelf if available else return OBSOLETE 
           and retrieve new block from queue  */
        next_block++;
        pop_record_by_block_id(records, next_block, &rec);
      } while (rec.block_id != RECORD_SLOT_OBSOLETE);
    } else {                    /* shelf the block if not next */
      put_into_record_v(records, rec);
    }
  }

  free_record_v(records);
  if (c->outfn) {    /* for stdout, will close at the end of main */
    fflush(out);
    fclose(out);
  }
  return 0;
}

#define episnp_test(snps, i) snps[(i)>>3]&(1<<((i)&0x7))
#define episnp_set(snps, i) snps[(i)>>3] |= 1<<((i)&0x7)

static void format_epiread(kstring_t *epi, bam1_t *b, refseq_t *rs, uint8_t bsstrand, char *chrm, window_t *w, uint8_t *snps, uint32_t snp_beg) {

  kstring_t es; int first_snp_loc = -1;
  es.s = 0; es.l = es.m = 0;
  
  int i; uint32_t j;
  bam1_core_t *c = &b->core;
  uint32_t rpos = c->pos+1, qpos = 0;
  int first_cpg_loc = -1;
  char qb, rb;
  for (i=0; i<c->n_cigar; ++i) {
    uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
    uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
    switch(op) {
    case BAM_CMATCH:
      for (j=0; j<oplen; ++j) {
        rb = toupper(getbase_refseq(rs, rpos+j));
        qb = bscall(b, qpos+j);
        if (bsstrand && rb == 'G' && rpos+j-1 >= rs->beg) {
          char rb0 = toupper(getbase_refseq(rs, rpos+j-1));
          if (rb0 == 'C') {	/* CpG context */
            if (first_cpg_loc < 0) {
              first_cpg_loc = (int) rpos+j-1;
              if ((unsigned) first_cpg_loc < w->beg || (unsigned) first_cpg_loc >= w->end)
                return;
              ksprintf(epi, "%s\t%d\t", chrm, first_cpg_loc);
            }
            if (qb == 'A') {
              kputc('T', epi);
            } else if (qb == 'G') {
              kputc('C', epi);
            } else {
              kputc('N', epi);
            }
          }
        }
        if (!bsstrand && rb == 'C' && rpos+j+1 <= rs->end) {
          char rb1 = toupper(getbase_refseq(rs, rpos+j+1));
          if (rb1 == 'G') {	/* CpG context */
            if (first_cpg_loc < 0) {
              first_cpg_loc = (int) rpos+j;
              if ((unsigned) first_cpg_loc < w->beg || (unsigned) first_cpg_loc >= w->end)
                return;
              ksprintf(epi, "%s\t%d\t", chrm, first_cpg_loc);
            }
            if (qb == 'T') {
              kputc('T', epi);
            } else if (qb == 'C') {
              kputc('C', epi);
            } else {
              kputc('N', epi);
            }
          }
        }

        /* append SNP info if present */
        uint32_t snp_ind = rpos+j-snp_beg;
        if (episnp_test(snps, snp_ind)) {
          kputc(qb, &es);
          if (first_snp_loc < 0)
            first_snp_loc = rpos+j;
        }
      }
      rpos += oplen;
      qpos += oplen;
      break;
    case BAM_CINS:
      qpos += oplen;
      break;
    case BAM_CDEL:
      rpos += oplen;
      break;
    case BAM_CSOFT_CLIP:
      qpos += oplen;
      break;
    case BAM_CHARD_CLIP:
      qpos += oplen;
      break;
    default:
      fprintf(stderr, "Unknown cigar, %u\n", op);
      abort();
    }
  }
  if (first_cpg_loc >= 0) {
    if (first_snp_loc >= 0)
      ksprintf(epi, "\t%d\t%s", first_snp_loc, es.s);
    kputc('\n', epi);
  }

  free(es.s);
}

static void *process_func(void *data) {

  result_t *res = (result_t*) data;
  conf_t *conf = (conf_t*) res->conf;
  samfile_t *in = samopen(res->bam_fn, "rb", 0);
  bam_index_t *idx = bam_index_load(res->bam_fn);
  refseq_t *rs = init_refseq(res->ref_fn, 1000, 1000);
  uint32_t j;

  record_t rec;
  memset(&rec, 0, sizeof(record_t));
  window_t w;
  while (1) {

    wqueue_get(window, res->q, &w);
    if (w.tid == -1) break;

    rec.tid = w.tid;
    char *chrm = in->header->target_name[w.tid];

    uint32_t snp_beg = w.beg>1000?w.beg-1000:1;
    uint32_t snp_end = w.end+1000;
    uint8_t *snps = calloc((snp_end-snp_beg)/8+1, sizeof(uint8_t)); /* TODO free snps */
    episnp_chrom1_t *episnp1 = get_episnp1(res->snp, chrm);
    for (j=0; j<episnp1->n; ++j) {
      uint32_t l=episnp1->locs[j];
      if (l>=snp_beg && l<snp_end) {
        episnp_set(snps, l-snp_beg);
      }
    }
    
    rec.s.l = rec.s.m = 0; rec.s.s = 0; /* the epiread string */

    fetch_refseq(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
    bam_iter_t iter = bam_iter_query(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
    bam1_t *b = bam_init1();
    int ret;
    while ((ret = bam_iter_read(in->x.bam, iter, b))>0) {

      uint8_t bsstrand = get_bsstrand(rs, b, conf->min_base_qual);

      /* read-based filtering */
      bam1_core_t *c = &b->core;
      if (c->qual < conf->min_mapq) continue;
      if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->min_read_len) continue;
      if (c->flag > 0){         /* only when any flag is set */
        if (conf->filter_secondary && c->flag & BAM_FSECONDARY) continue;
        if (conf->filter_duplicate && c->flag & BAM_FDUP) continue;
        if (conf->filter_ppair && !(c->flag & BAM_FPROPER_PAIR)) continue;
        if (conf->filter_qcfail && c->flag & BAM_FQCFAIL) continue;
      }

      uint8_t *nm = bam_aux_get(b, "NM");
      if (nm && bam_aux2i(nm)>conf->max_nm) continue;
      uint32_t cnt_ret = cnt_retention(rs, b, bsstrand);
      if (cnt_ret > conf->max_retention) continue;

      /* produce epiread */
      format_epiread(&rec.s, b, rs, bsstrand, chrm, &w, snps, snp_beg);
    }

    /* run through cytosines */
    rec.s.l = rec.s.m = 0; rec.s.s = 0; /* the record string */
    rec.block_id = w.block_id;

    /* put output string to output queue */
    wqueue_put2(record, res->rq, rec);

    bam_destroy1(b);
    bam_iter_destroy(iter);
  }
  free_refseq(rs);
  samclose(in);
  bam_index_destroy(idx);
  return 0;
}

static int usage(conf_t *conf) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: epiread [options] -r [ref.fa] -i [in.bam] -o [out.pileup] -g [chr1:123-234]\n");
  fprintf(stderr, "Input options:\n\n");
  fprintf(stderr, "     -i        input bam.\n");
  fprintf(stderr, "     -r        reference in fasta.\n");
  fprintf(stderr, "     -g        region (optional, if not specified the whole bam will be processed).\n");
  fprintf(stderr, "     -s        step of window dispatching [%d].\n", conf->step);
  fprintf(stderr, "     -q        number of threads [%d].\n", conf->n_threads);
  fprintf(stderr, "\nOutputing format:\n\n");
  fprintf(stderr, "     -o        output file [stdout]\n");
  fprintf(stderr, "     -N        NOMe-seq mode (skip GCH in bisulfite conversion estimate) [off]\n");
  fprintf(stderr, "     -R        epiread output file name [no]\n");
  fprintf(stderr, "     -B        bed input for SNP display in epiread output [no SNP]\n");
  fprintf(stderr, "     -v        verbose (print additional info for diagnosis).\n");
  fprintf(stderr, "\nPileup filtering:\n\n");
  fprintf(stderr, "     -k        min read coverage [%d]\n", conf->min_cov);
  fprintf(stderr, "     -b        min base quality [%u].\n", conf->min_base_qual);
  fprintf(stderr, "     -m        minimum mapping quality [%u].\n", conf->min_mapq);
  fprintf(stderr, "     -t        max cytosine retention in a read [%u].\n", conf->max_retention);
  fprintf(stderr, "     -l        minimum read length [%u].\n", conf->min_read_len);
  fprintf(stderr, "     -e        minimum distance to end of a read [%u].\n", conf->min_dist_end);
  fprintf(stderr, "     -c        NO filtering secondary mapping.\n");
  fprintf(stderr, "     -u        NO filtering of duplicate.\n");
  fprintf(stderr, "     -p        NO filtering of improper pair (!BAM_FPROPER_PAIR).\n");
  fprintf(stderr, "     -S        bsrate maximum position [%d]\n", conf->bsrate_max_pos);
  fprintf(stderr, "     -n        maximum NM tag [%u].\n", conf->max_nm);
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");

  return 1;
}

episnp_chrom1_v *bed_init_episnp(char *snp_bed_fn) {

  episnp_chrom1_v *episnp = init_episnp_chrom1_v(2);
  kstring_t line;
  line.l = line.m = 0; line.s = 0;

  episnp_chrom1_t *episnp1 = 0;
  char *tok;
  FILE *fh = fopen(snp_bed_fn,"r");
  while (1) {
    int c=fgetc(fh);
    if (c=='\n' || c==EOF) {
      tok = strtok(line.s, "\t");

      if (!episnp1 || strcmp(episnp1->chrm, tok) != 0)
        episnp1 = get_n_insert_episnp1(episnp, tok);

      episnp1->locs = realloc(episnp1->locs, (episnp1->n+1)*sizeof(uint32_t));
      tok = strtok(NULL, "\t");
      ensure_number(tok);
      episnp1->locs[episnp1->n] = atoi(tok)+1;
      episnp1->n++;
      
      line.l = 0;
      if (c==EOF) {
        break;
      }
    } else {
      kputc(c, &line);
    }
    free(line.s);
  }

  return episnp;
}

int main_epiread(int argc, char *argv[]) {

  int c;
  char *reffn = 0;
  char *reg = 0;
  char *infn = 0;
  char *outfn = 0;
  char *statsfn = 0;
  char *snp_bed_fn = 0;

  conf_t conf;
  memset(&conf, 0, sizeof(conf_t));
  conf.step = 100000;
  conf.n_threads = 3;
  conf.bsrate_max_pos = 1000;
  conf.min_base_qual = 20;
  conf.min_mapq = 40;
  conf.min_cov = 3;
  conf.max_retention = 999999;
  conf.min_read_len = 10;
  conf.filter_qcfail = 1;
  conf.filter_secondary = 1;
  conf.filter_duplicate = 1;
  conf.filter_ppair = 1;
  conf.min_dist_end = 3;
  conf.max_nm = 255;
  conf.is_nome = 0;
  conf.verbose = 0;

  if (argc<2) return usage(&conf);
  while ((c=getopt(argc, argv, "i:B:o:r:g:q:e:s:b:S:k:t:n:m:l:Ncupvh"))>=0) {
    switch (c) {
    case 'i': infn = optarg; break;
    case 'B': snp_bed_fn = optarg; break;
    case 'o': outfn = optarg; break;
    case 'r': reffn = optarg; break;
    case 'g': reg = optarg; break;
    case 'q': conf.n_threads = atoi(optarg); break;
    case 's': conf.step = atoi(optarg); break;
    case 'e': conf.min_dist_end = atoi(optarg); break;
    case 'b': conf.min_base_qual = atoi(optarg); break;
    case 'S': conf.bsrate_max_pos = atoi(optarg); break;
    case 'k': conf.min_cov = atoi(optarg); break;
    case 't': conf.max_retention = atoi(optarg); break;
    case 'l': conf.min_read_len = atoi(optarg); break;
    case 'n': conf.max_nm = atoi(optarg); break;
    case 'm': conf.min_mapq = atoi(optarg); break;
    case 'N': conf.is_nome = 1; break;
    case 'c': conf.filter_secondary = 0; break;
    case 'u': conf.filter_duplicate = 0; break;
    case 'p': conf.filter_ppair = 0; break;
    case 'v': conf.verbose = 1; break;
    case 'h': return usage(&conf);
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
  }

  if (!infn || !reffn) {
    usage(&conf);
    exit(1);
  }

  episnp_chrom1_v *episnp = NULL;
  if (conf.epiread && snp_bed_fn)
    episnp = bed_init_episnp(snp_bed_fn);

  wqueue_t(window) *wq = wqueue_init(window, 100000);
  pthread_t *processors = calloc(conf.n_threads, sizeof(pthread_t));
  result_t *results = calloc(conf.n_threads, sizeof(result_t));
  int i; unsigned j;
  samfile_t *in = samopen(infn, "rb", 0);

  /* sort sequence name by alphabetic order, chr1, chr10, chr11 ... */
  target_v *targets = init_target_v(50);
  target_t *t;
  for (i=0; i<in->header->n_targets; ++i) {
    t = next_ref_target_v(targets);
    t->tid = i;
    t->name = in->header->target_name[i];
    t->len = in->header->target_len[i];
  }
  qsort(targets->buffer, targets->size, sizeof(target_t), compare_targets);

  /* setup writer */
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
  for (i=0; i<conf.n_threads; ++i) {
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

  /* process bam */
  if (reg) {                    /* regional */
    int tid;
    uint32_t beg, end;
    bam_parse_region(in->header, reg, &tid, (int*) &beg, (int*) &end);
    /* chromosome are assumed to be less than 2**29 */
    beg++; end++;
    if (beg<=0) beg = 1;
    if (end>in->header->target_len[tid]) end = in->header->target_len[tid];
    for (wbeg = beg; wbeg < end; wbeg += conf.step, block_id++) {
      w.tid = tid;
      w.block_id = block_id;
      w.beg = wbeg;
      w.end = wbeg + conf.step;
      if (w.end > end) w.end = end;
      wqueue_put(window, wq, &w);
    }
  } else {                      /* entire bam */
    for (j=0; j<targets->size; ++j) {
      t = ref_target_v(targets, j);
      for (wbeg = 1; wbeg < t->len; wbeg += conf.step, block_id++) {
        w.tid = t->tid;
        w.block_id = block_id;
        w.beg = wbeg;
        w.end = wbeg+conf.step;
        if (w.end > t->len) w.end = t->len;
        wqueue_put(window, wq, &w);
      }
    }
  }
  for (i=0; i<conf.n_threads; ++i) {
    w.tid = -1;
    wqueue_put(window, wq, &w);
  }

  for (i=0; i<conf.n_threads; ++i) {
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
  samclose(in);

  return 0;
}
