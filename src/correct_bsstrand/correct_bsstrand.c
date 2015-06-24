/* 
 * correct bisulfite strand information if it is very inconsistent with C2T/G2A count
 * by zhouwanding@gmail.com
 *
 * if nC2T > 1 and nG2A > 1: then fail and mark "?"
 * if nC2T - nG2A > 2 and "-" => "+"
 * if nG2A - nC2T > 2 and "+" => "-"
 *
 * output a summary of how many reads are inconsistent
 *
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "sam.h"
#include "wstr.h"
#include "refseq.h"

typedef int (*sam_fetch_f)(bam1_t *b, const samfile_t *in, samfile_t *out, void *data);

int sam_fetch(char *ifn, char *ofn, char *reg, void *data, sam_fetch_f func) {
	int ret = 0;
	samfile_t *in = samopen(ifn, "rb", 0);
	samfile_t *out = 0;
	if (ofn) out = samopen(ofn, "wb", in->header);

	if (reg) {
		bam_index_t *idx = bam_index_load(ifn);
		if (idx == 0) {
			fprintf(stderr, "[%s:%d] Random alignment retrieval only works for indexed BAM files.\n",
							__func__, __LINE__);
			exit(1);
		}
		int tid, beg, end;
		bam_parse_region(in->header, reg, &tid, &beg, &end);
		if (tid < 0) {
			fprintf(stderr, "[%s:%d] Region \"%s\" specifies an unknown reference name. \n",
							__func__, __LINE__, reg);
			exit(1);
		}
		bam_iter_t iter;
		bam1_t *b = bam_init1();
		iter = bam_iter_query(idx, tid, beg, end);
		while ((ret = bam_iter_read(in->x.bam, iter, b)) >= 0) func(b, in, out, data);
		bam_iter_destroy(iter);
		bam_destroy1(b);
		bam_index_destroy(idx);
	} else {
		bam1_t *b = bam_init1();
		while ((ret = samread(in, b)) >= 0) func(b, in, out, data);
		bam_destroy1(b);
	}
	if (out) samclose(out);
	samclose(in);
			
	if (ret != -1) {					/* truncated is -2 */
		fprintf(stderr, "[%s:%d] Alignment retrieval failed due to truncated file\n",
						__func__, __LINE__);
		exit(1);
	}

	return ret;
}

#define bscall(seq, pos) bam_nt16_rev_table[bam1_seqi(seq, pos)]

typedef struct {
	uint8_t output_read, output_all_read;
	uint8_t infer_bsstrand;
} bsstrand_conf_t;

typedef struct {
	refseq_t *rs;
	int n_corr, n_mapped, n_fail, n_unmapped;
	bsstrand_conf_t *conf;
} bsstrand_data_t;

double similarity(int n1, int n2) {
	return ((double) (n1>n2? n2: n1) / (double) (n1>n2? n1 : n2));
}

int bsstrand_func(bam1_t *b, const samfile_t *in, samfile_t *out, void *data) {

	bsstrand_data_t *d = (bsstrand_data_t*)data;
	bsstrand_conf_t *conf = d->conf;
	const bam1_core_t *c = &b->core;

	if (c->flag & BAM_FUNMAP){
		if (out) samwrite(out, b);
		d->n_unmapped++;
		return 0;
	}
	
	fetch_refseq(d->rs, in->header->target_name[c->tid], c->pos, c->pos+1);
	uint32_t rpos=c->pos+1, qpos=0;
	int i, nC2T = 0, nG2A = 0;
	uint32_t j;
	char rbase, qbase;

	for (i=0; i<c->n_cigar; ++i) {
		uint32_t op = bam_cigar_op(bam1_cigar(b)[i]);
		uint32_t oplen = bam_cigar_oplen(bam1_cigar(b)[i]);
		switch(op) {
		case BAM_CMATCH:
			for(j=0; j<oplen; ++j) {
				rbase = toupper(getbase_refseq(d->rs, rpos+j));
				qbase = bscall(bam1_seq(b), qpos+j);
				if (rbase == 'C' && qbase == 'T') nC2T += 1;
				if (rbase == 'G' && qbase == 'A') nG2A += 1;
				/* printf("%c vs %c\n", toupper(rbase), qbase); */
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

	char key[2] = {'Z','S'};
	unsigned char *bsstrand = bam_aux_get(b, key);
	if (bsstrand) {
		bsstrand++;
		double s = similarity(nG2A, nC2T);
		if (nG2A > 1 && nC2T > 1 && s > 0.5) {
			if (conf->output_read || conf->output_all_read)
				printf("F\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", in->header->target_name[c->tid], c->pos, nC2T, nG2A, bam1_qname(b), bsstrand, s);
			bam_aux_append(b, "OS", 'A', 1, bsstrand);
			bsstrand[0] = '?';
			d->n_fail++;
		} else if (*bsstrand == '+' && nG2A > nC2T + 2) {
			if (conf->output_read || conf->output_all_read)
				printf("W2C\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", in->header->target_name[c->tid], c->pos, nC2T, nG2A, bam1_qname(b), bsstrand, s);
			bam_aux_append(b, "OS", 'A', 1, bsstrand);
			bsstrand[0] = '-';
			d->n_corr++;
		} else if (*bsstrand == '-' && nC2T > nG2A + 2) {
			if (conf->output_read || conf->output_all_read)
				printf("C2W\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", in->header->target_name[c->tid], c->pos, nC2T, nG2A, bam1_qname(b), bsstrand, s);
			bam_aux_append(b, "OS", 'A', 1, bsstrand);
			bsstrand[0] = '+';
			d->n_corr++;
		} else if (conf->output_all_read) {
			printf("N\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", in->header->target_name[c->tid], c->pos, nC2T, nG2A, bam1_qname(b), bsstrand, s);
		}
	} else if (!(c->flag & BAM_FUNMAP) && conf->infer_bsstrand) {
		char bss[3];
		if (similarity(nG2A, nC2T) < 0.5) {
			strcpy(bss, "??");
		} else if (nC2T > nG2A) {
			strcpy(bss, c->flag & BAM_FREVERSE ? "+-" : "++");
		} else {
			strcpy(bss, c->flag & BAM_FREVERSE ? "-+" : "--");
		}
		bam_aux_append(b, "ZS", 'Z', 3, (uint8_t*) bss);
	}

	
	if (out) samwrite(out, b);
	d->n_mapped++;

	return 0;
}

int bsstrand_main(char *reffn, char *ifn, char *ofn, bsstrand_conf_t *conf, char *reg) {

	bsstrand_data_t d = {
		.n_corr = 0,
		.n_mapped = 0,
		.n_fail = 0,
		.n_unmapped = 0,
		.rs = init_refseq(reffn, 100, 100000),
		.conf = conf,
	};
	sam_fetch(ifn, ofn, reg, &d, bsstrand_func);

	fprintf(stderr, "[%s:%d] Mapped reads: %d\n", __func__, __LINE__, d.n_mapped);
	fprintf(stderr, "[%s:%d] Unmapped reads: %d\n", __func__, __LINE__, d.n_unmapped);
	fprintf(stderr, "[%s:%d] Corrected reads: %d (%1.2f%%)\n",
					__func__, __LINE__, d.n_corr, (double)d.n_corr/(double)d.n_mapped*100.);
	fprintf(stderr, "[%s:%d] Failed reads: %d (%1.2f%%)\n",
					__func__, __LINE__, d.n_fail, (double)d.n_fail/(double)d.n_mapped*100.);
	free_refseq(d.rs);
	return 0;
}


static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: correct_bsstrand [options] -r [ref.fa] -i [in.bam] -o [out.bam] -e [readlist] -g [chr1:123-234]\n");
  fprintf(stderr, "Input options:\n");
	fprintf(stderr, "     -i        input bam.\n");
  fprintf(stderr, "     -r        reference in fasta.\n");
	fprintf(stderr, "     -g        region (optional, if not specified the whole bam will be processed).\n");
	fprintf(stderr, "     -o        output bam (optional).\n");
	fprintf(stderr, "     -e        output abnormal read stats [optional].\n");
	fprintf(stderr, "     -a        output all reads stats [optional].\n");
	fprintf(stderr, "     -m        infer bsstrand information if missing [optional].\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {
	int c;
	char *reffn = 0;
	char *reg = 0;										/* region */
	char *infn = 0;
	char *outfn = 0;
	bsstrand_conf_t conf = {.output_read = 0, .output_all_read = 0, .infer_bsstrand = 0};

  if (argc < 2) return usage();
  while ((c = getopt(argc, argv, "i:o:r:g:eamh")) >= 0) {
    switch (c) {
		case 'i': infn = optarg; break;
		case 'o': outfn = optarg; break;
		case 'r': reffn = optarg; break;
		case 'g': reg = optarg; break;
    case 'e': conf.output_read = 1; break;
		case 'a': conf.output_all_read = 1; break;
		case 'm': conf.infer_bsstrand = 1; break;
    case 'h': return usage();
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
  }

	bsstrand_main(reffn, infn, outfn, &conf, reg);

	return 0;
}
