#include <ctype.h>
#include "wstr.h"
#include "wqueue.h"
#include "encode.h"
#include "sam.h"
#include "refseq.h"
#include "kstring.h"

typedef struct {
	int32_t tid;
  uint32_t beg, end;
} window_t;

DEFINE_WQUEUE(window, window_t)

uint8_t trinuc2ind(char *trinuc) {
	uint8_t c = 0;
	/* if (trinuc[1] == 'C') c = 0; */
	/* else if (trinuc[1] == 'G') c = 1; */
	/* else abort(); */
	c = (c<<2)|nt256char_to_nt256int8_table[(unsigned char)trinuc[0]];
	c = (c<<2)|nt256char_to_nt256int8_table[(unsigned char)trinuc[2]];
	return c;
}

void ind2trinuc(char *trinuc, uint8_t c) {
	trinuc[0] = nt256int8_to_nt256char_table[(c>>2)&3];
	trinuc[2] = nt256int8_to_nt256char_table[c&3];
	trinuc[1] = 'C';
	/* trinuc[1] = c&(1>>4)?'G':'C'; */
}

typedef struct {
	int step;
	int n_threads;
	char *output_retention;
	uint32_t min_base_qual;
	uint32_t max_retention;
	uint8_t min_dist_end;
	uint8_t relax;
} conf_t;

DEFINE_WQUEUE(record, char*)

/* typedef struct { */
/* 	bam1_t *b; */
/* 	int rpos; */
/* 	int qpos; */
/* 	int bsstrand; */
/* 	char context[3]; */
/* } record_t; */


/* DEFINE_WQUEUE(record, record_t) */

typedef struct {
  int n_targets;
	char *bam_fn;									/* on stack */
	char *ref_fn;									/* on stack */
  long **total;
  long **retained;
  wqueue_t(window) *q;
	wqueue_t(record) *rq;
	conf_t *conf;
} result_t;

result_t *init_result(int n_targets) {

	int i;
  result_t *result = (result_t*) calloc(1,sizeof(result_t));
  result->n_targets = n_targets;
  result->total = (long**) malloc(n_targets*sizeof(long*));
  result->retained = (long**) malloc(n_targets*sizeof(long*));
  for (i=0; i<n_targets; ++i) {
    result->total[i] = (long*) calloc(16, sizeof(long));
    result->retained[i] = (long*) calloc(16, sizeof(long));
  }
  return result;
}

typedef struct {
	wqueue_t(record) *q;
	char *outfn;
} writer_conf_t;

/* void *write_func(void *data) { */
/* 	writer_conf_t *c = (writer_conf_t*) data; */
/* 	FILE *out = fopen(c->outfn, "w"); */
/* 	while (1) { */
/* 		record_t rec; */
/* 		wqueue_get(record, c->q, &rec); */
/* 		if(!rec.b) break; */
/* 		fprintf(out, "%s\t%d\t%d\t%c%c%c\n", bam1_qname(rec.b), rec.rpos, rec.qpos, rec.context[0], rec.context[1], rec.context[2]); */
/* 		bam_destroy1(rec.b); */
/* 	} */
/* 	fclose(out); */
/* 	return 0; */
/* } */

void *write_func(void *data) {
	writer_conf_t *c = (writer_conf_t*) data;
	FILE *out = fopen(c->outfn, "w");
	fputs("chrm\tpos\tbsstrand\tcontext\tstrand\tpos_on_reads\tqname\tread_len\n", out);
	while (1) {
		char *rec;
		wqueue_get(record, c->q, &rec);
		if(!rec) break;
		fputs(rec, out);
		free(rec);									/* rec is alloc-ed */
	}
	fclose(out);
	return 0;
}

void destroy_result(result_t *result) {
	int i;
  for (i=0; i<result->n_targets; ++i) {
    free(result->total[i]);
    free(result->retained[i]);
  }
  free(result->total);
  free(result->retained);
  free(result);
}

void merge_results(result_t *dst, result_t *src) {
	int i, j;
	for (i=0; i<src->n_targets; ++i) {
		for (j=0; j<16; ++j) {
			dst->total[i][j] += src->total[i][j];
			dst->retained[i][j] += src->retained[i][j];
		}
	}
}

#define bscall(seq, pos) bam_nt16_rev_table[bam1_seqi(seq, pos)]

void profile_read(bam1_t *b, refseq_t *rs, uint8_t *bsstrand, uint32_t *cnt_retention, uint32_t *nC2T, uint32_t *nG2A) {

	*cnt_retention = 0; *nC2T = 0; *nG2A = 0;
	int i; uint32_t j;
	bam1_core_t *c = &b->core;
	uint32_t rpos =c->pos+1, qpos = 0;
	for (i=0; i<c->n_cigar; ++i) {
		uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
		uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
		char rb, qb;
		switch(op) {
		case BAM_CMATCH:
			for (j=0; j<oplen; ++j) {
				rb = getbase_refseq(rs, rpos+j);
				qb = bscall(bam1_seq(b), qpos+j);
				if (rb == 'C' && qb == 'T') (*nC2T)++;
				if (rb == 'G' && qb == 'A') (*nG2A)++;
				if (bsstrand[0] == '+' && rb == 'C' && qb == 'C') (*cnt_retention)++;
				if (bsstrand[0] == '-' && rb == 'G' && qb == 'G') (*cnt_retention)++;
				else  continue;
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
		default:
			fprintf(stderr, "Unknown cigar, %u\n", op);
			abort();
		}
	}
}


void *process_func(void *data) {

  result_t *res = (result_t*) data;
	conf_t *conf = (conf_t*) res->conf;
	samfile_t *in = samopen(res->bam_fn, "rb", 0);
	bam_index_t *idx = bam_index_load(res->bam_fn);
	refseq_t *rs = init_refseq(res->ref_fn, 1000, 1000);

	char qbase;
	char *rbase;
	char trinuc[3];
	char trinuc_r[3];
	window_t w;

	int i; uint32_t j;
  while (1) {

    wqueue_get(window, res->q, &w);
    if (w.tid == -1) break;
		/* since we need to look ahead and after */
		char *chrm = in->header->target_name[w.tid];
		fetch_refseq(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
		bam_iter_t iter = bam_iter_query(idx, w.tid, w.beg, w.end);
		bam1_t *b = bam_init1();
		int ret;
		while ((ret = bam_iter_read(in->x.bam, iter, b))>0) {

			uint8_t *bsstrand = bam_aux_get(b, "ZS");
			if (!bsstrand) continue;
			bsstrand++;

			bam1_core_t *c = &b->core;
			uint32_t rpos = c->pos+1, qpos = 0;
			if (rpos <= w.beg) continue;
			if (c->flag & BAM_FSECONDARY) continue;
			uint8_t *nm = bam_aux_get(b, "NM");
			if (nm && bam_aux2i(nm)>2) continue;

			uint32_t cnt_retention = 0, nC2T = 0, nG2A = 0;
			profile_read(b, rs, bsstrand, &cnt_retention, &nC2T, &nG2A);
			if (!conf->relax) {
				if (nC2T > 0 && nG2A > 0) continue;
				if (nC2T > 0 && bsstrand[0] == '-') bsstrand[0] = '+';
				if (nG2A > 0 && bsstrand[0] == '+') bsstrand[0] = '-';
			}
			if (cnt_retention >= conf->max_retention) continue;

			for (i=0; i<c->n_cigar; ++i) {
				uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
				uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
				switch(op) {
				case BAM_CMATCH:
					for (j=0; j<oplen; ++j) {
						if(qpos+j <= conf->min_dist_end ||
							 c->l_qseq-qpos-j <= conf->min_dist_end) continue;
						if(bam1_qual(b)[qpos+j] <= conf->min_base_qual) continue;
						rbase = subseq_refseq(rs, rpos+j);
						if (rpos+j <= 1) continue; /* skip first base */
						if (rpos+j >= in->header->target_len[w.tid]) continue; /* skip last base */
						trinuc[0] = toupper(*(rbase-1));
						trinuc[1] = toupper(*rbase);
						trinuc[2] = toupper(rbase[1]);
						if (trinuc[0] == 'N' || trinuc[1] == 'N' || trinuc[2] == 'N')
							continue;

						qbase = bscall(bam1_seq(b), qpos+j);
						// printf("%c\t%c\n", trinuc[1], qbase);
						if (bsstrand[0] == '+' && trinuc[1] == 'C') {
							res->total[c->tid][trinuc2ind(trinuc)]++;
							if (qbase == 'C') {
								res->retained[c->tid][trinuc2ind(trinuc)]++;
								if (res->rq) {
									wqueue_put2(record, res->rq,
															wasprintf("%s\t%u\t%c\t%.3s\t%c\t%u\t%s\t%d\n",
																				chrm, rpos, bsstrand[0], trinuc,
																				c->flag & BAM_FREVERSE ? '-' : '+',
																				c->flag & BAM_FREVERSE ? c->l_qseq-qpos-j : qpos+j+1,
																				bam1_qname(b), c->l_qseq));
								}
							}
						} else if (bsstrand[0] == '-' && trinuc[1] == 'G') {
							_nt256char_rev(trinuc_r, trinuc, 3);
							res->total[c->tid][trinuc2ind(trinuc_r)]++;
							if (qbase == 'G') {
								res->retained[c->tid][trinuc2ind(trinuc_r)]++;
								if (res->rq) {
									wqueue_put2(record, res->rq,
															wasprintf("%s\t%u\t%c\t%.3s\t%c\t%u\t%s\t%d\n",
																				chrm, rpos, bsstrand[0], trinuc,
																				c->flag & BAM_FREVERSE ? '-' : '+',
																				c->flag & BAM_FREVERSE ? c->l_qseq-qpos-j : qpos+j+1,
																				bam1_qname(b), c->l_qseq));
								}
							}
						} else {
							continue;
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
				default:
					fprintf(stderr, "Unknown cigar, %u\n", op);
					abort();
				}
			}
		}
		bam_destroy1(b);
		bam_iter_destroy(iter);
  }
	free_refseq(rs);
	samclose(in);
	bam_index_destroy(idx);
	return 0;
}


int sample_trinuc_main(char *reffn, char *infn, char *reg, conf_t *conf) {

  wqueue_t(window) *wq = wqueue_init(window, 100000);
  pthread_t *processor = (pthread_t*) calloc(conf->n_threads, sizeof(pthread_t));
  result_t **results = (result_t**) malloc(conf->n_threads*sizeof(result_t*));
	int i, j;
	samfile_t *in = samopen(infn, "rb", 0);
	pthread_t writer;
	writer_conf_t writer_conf = {.q = 0, .outfn = 0};
	if (conf->output_retention) {
		writer_conf.q = wqueue_init(record, 100000);
		writer_conf.outfn = conf->output_retention;
		pthread_create(&writer, NULL, write_func, &writer_conf);
	}
  for (i=0; i<conf->n_threads; ++i) {
    results[i] = init_result(in->header->n_targets);
		results[i]->q = wq;
		results[i]->rq = writer_conf.q;
		results[i]->ref_fn = reffn;
		results[i]->bam_fn = infn;
		results[i]->conf = conf;
    pthread_create(&processor[i], NULL, process_func, results[i]);
  }

	window_t w; memset(&w, 0, sizeof(window_t));
	uint32_t wbeg;
	if (reg) {
		int tid;
		uint32_t beg, end;
		bam_parse_region(in->header, reg, &tid, (int*) &beg, (int*) &end);
		/* chromosome are assumed to be less than 2**29 */
		if (end>in->header->target_len[tid]) end = in->header->target_len[tid];
		for (wbeg = beg; wbeg < end; wbeg += conf->step) {
			w.tid = tid;
			w.beg = wbeg;
			w.end = wbeg + conf->step;
			if (w.end > end) w.end = end;
			wqueue_put(window, wq, &w);
		}
	} else {
		for (i=0; i<in->header->n_targets; ++i) {
			uint32_t target_len = in->header->target_len[i];
			for (wbeg = 1; wbeg < target_len; wbeg += conf->step) {
				w.tid = i;
				w.beg = wbeg;
				w.end = wbeg+conf->step;
				if (w.end > target_len) w.end = target_len;
				wqueue_put(window, wq, &w);
			}
		}
	}
  for (i=0; i<conf->n_threads; ++i) {
    w.tid = -1;
    wqueue_put(window, wq, &w);
  }

  for (i=0; i<conf->n_threads; ++i) {
    pthread_join(processor[i], NULL);
	}

	if (conf->output_retention) {
		wqueue_put2(record, writer_conf.q, NULL);
		/* record_t rec = {.b=0}; */
		/* wqueue_put(record, writer_conf.q, &rec); */
		pthread_join(writer, NULL);
		wqueue_destroy(record, writer_conf.q);
	}

  result_t *final_result = init_result(in->header->n_targets);
  for (i=0; i<conf->n_threads; ++i) {
    merge_results(final_result, results[i]);
    destroy_result(results[i]);
  }
  free(results);

  kstring_t s;
  s.l = s.m = 0; s.s = 0;
  kputs("trinuc", &s);
	char trinuc[3];
  for (j=0; j<16; ++j) {
    kputs("\tn", &s);
		ind2trinuc(trinuc, j);
    kputsn(trinuc, 3, &s);
  }
	for(j=0; j<16; ++j) {
		kputs("\tr", &s);
		ind2trinuc(trinuc,j);
		kputsn(trinuc, 3, &s);
	}
  puts(s.s);
  for (i=0; i<in->header->n_targets; ++i) {
		long total_all=0;
    s.l = s.m = 0; free(s.s); s.s = 0;
    kputs(in->header->target_name[i], &s);
    for (j=0; j<16; ++j) {
			ksprintf(&s, "\t%d", final_result->total[i][j]);
			total_all += final_result->total[i][j];
		}
    for (j=0; j<16; ++j) ksprintf(&s, "\t%d", final_result->retained[i][j]);
    if(total_all) puts(s.s);
  }
	if (s.s) free(s.s);
	destroy_result(final_result);

	free(processor);
  wqueue_destroy(window, wq);
	samclose(in);

	return 0;
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: sample_trinuc [options] -r [ref.fa] -i [in.bam] -g [chr1:123-234]\n");
  fprintf(stderr, "Input options:\n");
	fprintf(stderr, "     -i        input bam.\n");
  fprintf(stderr, "     -r        reference in fasta.\n");
	fprintf(stderr, "     -g        region (optional, if not specified the whole bam will be processed).\n");
	fprintf(stderr, "     -e        output retained read stats [false]\n");
	fprintf(stderr, "     -s        step of window dispatching [100000].\n");
	fprintf(stderr, "     -q        number of threads [3].\n");
	fprintf(stderr, "     -b        min base quality [10].\n");
	fprintf(stderr, "     -d        minimum distance to end [2].\n");
	fprintf(stderr, "     -t        max retention in a read [999999].\n");
	fprintf(stderr, "     -z        drop stringency filter, allow reads with nG2A>0 and nC2T>0 and no correction of bsstrand if nG2A/nC2T is not consistent with bsstrand.\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {

	int c;
	char *reffn = 0;
	char *reg = 0;
	char *infn = 0;
	conf_t conf = {
		.step = 100000,
		.n_threads = 3,
		.output_retention = NULL,
		.min_base_qual = 10,
		.max_retention = 999999,
		.min_dist_end = 2,
		.relax = 0,
	};

	if (argc<2) return usage();
	while ((c=getopt(argc, argv, "i:o:r:g:s:q:e:b:t:d:z:h"))>=0) {
		switch (c) {
		case 'i': infn = optarg; break;
		case 'r': reffn = optarg; break;
		case 'g': reg = optarg; break;
		case 's': conf.step = atoi(optarg); break;
		case 'e': conf.output_retention = optarg; break;
		case 'q': conf.n_threads = atoi(optarg); break;
		case 'b': conf.min_base_qual = atoi(optarg); break;
		case 't': conf.max_retention = atoi(optarg); break;
		case 'd': conf.min_dist_end = atoi(optarg); break;
		case 'z': conf.relax = 1; break;
		case 'h': return usage();
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
	}

	if (!infn || !reffn) usage();

	sample_trinuc_main(reffn, infn, reg, &conf);
}

