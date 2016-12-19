#include <ctype.h>
#include "wqueue.h"
#include "encode.h"
#include "sam.h"
#include "refseq.h"
#include "kstring.h"
#include "wvec.h"
#include "wstr.h"

typedef struct {
	int32_t tid;
  uint32_t beg, end;
} window_t;

DEFINE_WQUEUE(window, window_t)

typedef struct {
	int step;
	int n_threads;
	uint32_t min_base_qual;
	uint32_t max_retention;
	uint32_t min_read_len;
	uint8_t min_dist_end;
	uint32_t max_plpsize;
	uint8_t max_nm;
	uint8_t filter_secondary:1;
	uint8_t print_retention_only:1;
} conf_t;

typedef enum {BSS_RETENTION, BSS_CONVERSION, BSS_OTHER} bsstate_t;

typedef struct {
	char *qname;
	uint32_t beg;
	uint32_t end;
	char bsstrand;
	bsstate_t bsstate;
} insert_t;

DEFINE_VECTOR(insert_v, insert_t)

DEFINE_WQUEUE(record, char*)

typedef struct {
	char *bam_fn;									/* on stack */
	char *ref_fn;									/* on stack */
  wqueue_t(window) *q;
	wqueue_t(record) *rq;
	conf_t *conf;
} result_t;

typedef struct {
	wqueue_t(record) *q;
	char *outfn;
} writer_conf_t;

void *write_func(void *data) {
	writer_conf_t *c = (writer_conf_t*) data;
	FILE *out = fopen(c->outfn, "w");
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

#define bscall(b, pos) bam_nt16_rev_table[bam1_seqi(bam1_seq(b), pos)]

/* return -1 if abnormal (missing bsstrand) */
int cnt_retention(refseq_t *rs, bam1_t *b) {
	int cnt = 0;
	uint8_t *bsstrand = bam_aux_get(b, "ZS");
	if (!bsstrand) return -1;
	bsstrand++;

	bam1_core_t *c = &b->core;
	uint32_t rpos = c->pos+1, qpos = 0;
	uint32_t op, oplen;
	char rb, qb;
	int i; unsigned j;
	for (i=0; i<c->n_cigar; ++i) {
		op = bam_cigar_op(bam1_cigar(b)[i]);
		oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
		switch(op) {
		case BAM_CMATCH:
			for (j=0; j<oplen; ++j) {
				rb = toupper(getbase_refseq(rs, rpos+j));
				qb = bscall(b, qpos+j);
				if (bsstrand[0] == '+' && rb == 'C' && qb == 'C') cnt++;
				else if (bsstrand[0] == '-' && rb == 'G' && qb == 'G') cnt++;
				else continue;
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

	return cnt;
}

void *process_func(void *data) {

  result_t *res = (result_t*) data;
	conf_t *conf = (conf_t*) res->conf;
	samfile_t *in = samopen(res->bam_fn, "rb", 0);
	bam_index_t *idx = bam_index_load(res->bam_fn);
	refseq_t *rs = init_refseq(res->ref_fn, 1000, 1000);

	window_t w;
	int i; uint32_t j, k;
  while (1) {

    wqueue_get(window, res->q, &w);
    if (w.tid == -1) break;

		int n = w.end - w.beg;
		insert_v **pos2inserts = calloc(n, sizeof(insert_v*));
		
		char *chrm = in->header->target_name[w.tid];
		fetch_refseq(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
		bam_iter_t iter = bam_iter_query(idx, w.tid, w.beg, w.end);
		bam1_t *b = bam_init1();
		int ret;
		char qb, rb, rb_next, rb_prev;
		while ((ret = bam_iter_read(in->x.bam, iter, b))>0) {

			uint8_t *bsstrand = bam_aux_get(b, "ZS");
			if (!bsstrand) continue;
			bsstrand++;

			bam1_core_t *c = &b->core;

			if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->min_read_len) continue;

			uint32_t rpos = c->pos+1, qpos = 0;
			if (conf->filter_secondary && c->flag & BAM_FSECONDARY) continue;
			if (!(c->flag & BAM_FPROPER_PAIR)) continue;
			uint8_t *nm = bam_aux_get(b, "NM");
			if (nm && bam_aux2i(nm)>conf->max_nm) continue;

			int cnt_ret = cnt_retention(rs, b);
			if (cnt_ret < 0 || (unsigned) cnt_ret > conf->max_retention) continue;

			rpos = c->pos+1; qpos = 0;
			for (i=0; i<c->n_cigar; ++i) {
				uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
				uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
				switch(op) {
				case BAM_CMATCH:
					for (j=0; j<oplen; ++j) {
						rb = toupper(getbase_refseq(rs, rpos+j));
						qb = bscall(b, qpos+j);
						if (rpos+j<w.beg || rpos+j>w.end) continue; /* include begin but not end */
						if (rpos+j == w.end && qb != 'G') continue; /* cross window C|G */
						if (rb == 'C' && bsstrand[0] == '+') {
							if (rpos+j+1 >= (unsigned) rs->seqlen) continue;
							rb_next = toupper(getbase_refseq(rs, rpos+j+1));
							if (rb_next != 'G') continue;
							insert_v **inses = pos2inserts + rpos + j - w.beg;
							if (!*inses) *inses = init_insert_v(2);
							insert_t *insert = next_ref_insert_v(*inses);
							insert->qname = strdup(bam1_qname(b));
							insert->beg = c->pos<c->mpos?c->pos:c->mpos+1;
							insert->end = insert->beg + abs(c->isize);
							insert->bsstrand = bsstrand[0];
							if (qb == 'C') insert->bsstate = BSS_RETENTION;
							else if (qb == 'T') insert->bsstate = BSS_CONVERSION;
							else insert->bsstate = BSS_OTHER;
						} else if (rb == 'G' && bsstrand[0] == '-') {
							if (rpos+j-1 <= w.beg) continue;
							rb_prev = toupper(getbase_refseq(rs, rpos+j-1));
							if (rb_prev != 'C') continue;
							insert_v **inses = pos2inserts + rpos + j - w.beg - 1;
							if (!*inses) *inses = init_insert_v(2);
							insert_t *insert = next_ref_insert_v(*inses);
							insert->qname = strdup(bam1_qname(b));
							insert->beg = c->pos<c->mpos?c->pos:c->mpos+1;
							insert->end = insert->beg + abs(c->isize);
							insert->bsstrand = bsstrand[0];
							if (qb == 'G') insert->bsstate = BSS_RETENTION;
							else if (qb == 'A') insert->bsstate = BSS_CONVERSION;
							else insert->bsstate = BSS_OTHER;
						} else continue;
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

		for (j=w.beg; j<w.end; ++j) {
			rb = getbase_refseq(rs, j);
			insert_v *inses = pos2inserts[j-w.beg];
			if (inses && inses->size <= conf->max_plpsize) {
				unsigned i_bsw, i_bsc;
				for (i_bsw = 0; i_bsw < inses->size; ++i_bsw) {

					insert_t *bsw = ref_insert_v(inses, i_bsw);
					if (bsw->bsstrand != '+') continue;
					for (i_bsc = 0; i_bsc < inses->size; ++i_bsc) {
						insert_t *bsc = ref_insert_v(inses, i_bsc);
						if (bsc->bsstrand != '-') continue;
						if (abs(bsw->beg - bsc->beg) + abs(bsw->end - bsc->end)>1) continue;
						if ((bsw->bsstate == BSS_CONVERSION && bsc->bsstate == BSS_RETENTION) ||
								(bsw->bsstate == BSS_RETENTION && bsc->bsstate == BSS_CONVERSION)) {
							wqueue_put2(record, res->rq,
													wasprintf("%s\t%u\t%s\t%u-%u\t%d\t%s\t%u-%u\t%d\n",
																		chrm, j, bsw->qname, bsw->beg, bsw->end, bsw->bsstate, 
																		bsc->qname, bsc->beg, bsc->end, bsc->bsstate));
						}
					}
				}
			}
		}
		for (i=0; i<n; ++i) {
			if (pos2inserts[i]) {
				for (k=0; k<pos2inserts[i]->size; ++k) {
					insert_t *ins = ref_insert_v(pos2inserts[i], k);
					free(ins->qname);
				}
				free_insert_v(pos2inserts[i]);
			}
		}
		free(pos2inserts);
		bam_destroy1(b);
		bam_iter_destroy(iter);
  }
	free_refseq(rs);
	samclose(in);
	bam_index_destroy(idx);
	return 0;
}


int hemifinder_main(char *reffn, char *infn, char *outfn, char *reg, conf_t *conf) {

  wqueue_t(window) *wq = wqueue_init(window, 100000);
  pthread_t *processors = calloc(conf->n_threads, sizeof(pthread_t));
  result_t *results = calloc(conf->n_threads, sizeof(result_t));
	int i;
	samfile_t *in = samopen(infn, "rb", 0);

	pthread_t writer;
	writer_conf_t writer_conf = {
		.q = wqueue_init(record, 100000),
		.outfn = outfn,
	};
	pthread_create(&writer, NULL, write_func, &writer_conf);
  for (i=0; i<conf->n_threads; ++i) {
		results[i].q = wq;
		results[i].rq = writer_conf.q;
		results[i].ref_fn = reffn;
		results[i].bam_fn = infn;
		results[i].conf = conf;
    pthread_create(&processors[i], NULL, process_func, &results[i]);
  }

	window_t w; memset(&w, 0, sizeof(window_t));
	uint32_t wbeg;
	if (reg) {
		int tid;
		uint32_t beg, end;
		bam_parse_region(in->header, reg, &tid, (int*) &beg, (int*) &end);
		/* chromosome are assumed to be less than 2**29 */
		beg++; end++;
		if (beg<=0) beg = 1;
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
    pthread_join(processors[i], NULL);
	}

	wqueue_put2(record, writer_conf.q, NULL);
	pthread_join(writer, NULL);
	wqueue_destroy(record, writer_conf.q);

	free(results);
	free(processors);
  wqueue_destroy(window, wq);
	samclose(in);

	return 0;
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: hemifinder [options] -r [ref.fa] -i [in.bam] -o [out.pileup] -g [chr1:123-234]\n");
	fprintf(stderr, "output format: chrm, pos, pos, refbase, context, coverage, filtered_retention, filtered_conversion, num_retention_in_reads, strands, position_on_reads, quality, bsstrand\n");
  fprintf(stderr, "Input options:\n");
	fprintf(stderr, "     -i        input bam.\n");
  fprintf(stderr, "     -r        reference in fasta.\n");
	fprintf(stderr, "     -g        region (optional, if not specified the whole bam will be processed).\n");
	fprintf(stderr, "     -o        pileup output file\n");
	fprintf(stderr, "     -s        step of window dispatching [100000].\n");
	fprintf(stderr, "     -q        number of threads [3] recommend 20.\n");
	fprintf(stderr, "     -b        min base quality [10].\n");
	fprintf(stderr, "     -t        max retention in a read [999999].\n");
	fprintf(stderr, "     -l        minimum read length [10].\n");
	fprintf(stderr, "     -e        minimum distance to end of a read [3].\n");
	fprintf(stderr, "     -c        turn off filtering secondary mapping.\n");
	fprintf(stderr, "     -m        maximum size of pileup size [1000].\n");
	fprintf(stderr, "     -n        maximum nm tag [2].\n");
	fprintf(stderr, "     -a        print all, not just retention reads.\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {

	int c;
	char *reffn = 0;
	char *reg = 0;
	char *infn = 0;
	char *outfn = 0;
	conf_t conf = {
		.step = 100000,
		.n_threads = 3,
		.min_base_qual = 10,
		.max_retention = 999999,
		.min_read_len = 10,
		.filter_secondary = 1,
		.print_retention_only = 1,
		.max_plpsize = 1000,
		.min_dist_end = 3,
		.max_nm = 2,
	};


	if (argc<2) return usage();
	while ((c=getopt(argc, argv, "i:o:r:g:q:e:b:t:m:l:cah"))>=0) {
		switch (c) {
		case 'i': infn = optarg; break;
		case 'r': reffn = optarg; break;
		case 'g': reg = optarg; break;
		case 'o': outfn = optarg; break;
		case 's': conf.step = atoi(optarg); break;
		case 'q': conf.n_threads = atoi(optarg); break;
		case 'b': conf.min_base_qual = atoi(optarg); break;
		case 't': conf.max_retention = atoi(optarg); break;
		case 'l': conf.min_read_len = atoi(optarg); break;
		case 'e': conf.min_dist_end = atoi(optarg); break;
		case 'c': conf.filter_secondary = 0; break;
		case 'm': conf.max_plpsize = atoi(optarg); break;
		case 'a': conf.print_retention_only = 0; break;
		case 'n': conf.max_nm = atoi(optarg); break;
		case 'h': return usage();
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
	}

	if (!infn || !reffn || !outfn) usage();

	hemifinder_main(reffn, infn, outfn, reg, &conf);
}

