/* 
 * Correct bisulfite strand information if it is very inconsistent with C2T/G2A count
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
 *
 */

#include <unistd.h>
#include "wstr.h"
#include "wzmisc.h"
#include "refcache.h"
#include "sam.h"
#include "bamfilter.h"

typedef struct {
	uint8_t output_read, output_all_read;
	uint8_t infer_bsstrand;
} bsstrand_conf_t;

typedef struct {
	refcache_t *rs;
	int n_corr, n_mapped, n_fail, n_unmapped;
	bsstrand_conf_t *conf;
} bsstrand_data_t;

int bsstrand_func(bam1_t *b, const samFile *in, samFile *out, bam_hdr_t *header, void *data) {

	bsstrand_data_t *d = (bsstrand_data_t*)data;
	bsstrand_conf_t *conf = d->conf;
	const bam1_core_t *c = &b->core;

	if (c->flag & BAM_FUNMAP){
		if (out) 
      if (sam_write1(out, header, b) < 0)
        wzfatal("Cannot write bam.\n");
		d->n_unmapped++;
		return 0;
	}
	
	fetch_refcache(d->rs, header->target_name[c->tid], c->pos, bam_endpos(b)+1);
	uint32_t rpos=c->pos+1, qpos=0;
	int i, nC2T = 0, nG2A = 0;
	uint32_t j;
	char rbase, qbase;

	for (i=0; i<c->n_cigar; ++i) {
		uint32_t op = bam_cigar_op(bam_get_cigar(b)[i]);
		uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
		switch(op) {
		case BAM_CMATCH:
			for(j=0; j<oplen; ++j) {
				rbase = toupper(getbase_refcache(d->rs, rpos+j));
				qbase = bscall(b, qpos+j);
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
    case BAM_CHARD_CLIP:
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
		double s = min(nG2A, nC2T)/max(nG2A, nC2T);
		if (nG2A > 1 && nC2T > 1 && s > 0.5) {
			if (conf->output_read || conf->output_all_read)
				printf("F\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", header->target_name[c->tid], c->pos, nC2T, nG2A, bam_get_qname(b), bsstrand, s);
			bam_aux_append(b, "OS", 'A', 1, bsstrand);
			bsstrand[0] = '?';
			d->n_fail++;
		} else if (*bsstrand == '+' && nG2A > nC2T + 2) {
			if (conf->output_read || conf->output_all_read)
				printf("W2C\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", header->target_name[c->tid], c->pos, nC2T, nG2A, bam_get_qname(b), bsstrand, s);
			bam_aux_append(b, "OS", 'A', 1, bsstrand);
			bsstrand[0] = '-';
			d->n_corr++;
		} else if (*bsstrand == '-' && nC2T > nG2A + 2) {
			if (conf->output_read || conf->output_all_read)
				printf("C2W\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", header->target_name[c->tid], c->pos, nC2T, nG2A, bam_get_qname(b), bsstrand, s);
			bam_aux_append(b, "OS", 'A', 1, bsstrand);
			bsstrand[0] = '+';
			d->n_corr++;
		} else if (conf->output_all_read) {
			printf("N\t%s\t%d\t%d\t%d\t%s\t%s\t%1.2f\n", header->target_name[c->tid], c->pos, nC2T, nG2A, bam_get_qname(b), bsstrand, s);
		}
	} else if (!(c->flag & BAM_FUNMAP) && conf->infer_bsstrand) {
		char bss[3];
		if (min(nG2A, nC2T) / max(nG2A, nC2T) < 0.5) {
			strcpy(bss, "??");
		} else if (nC2T > nG2A) {
			strcpy(bss, c->flag & BAM_FREVERSE ? "+-" : "++");
		} else {
			strcpy(bss, c->flag & BAM_FREVERSE ? "-+" : "--");
		}
		bam_aux_append(b, "ZS", 'Z', 3, (uint8_t*) bss);
	}

	if (out) 
    if (sam_write1(out, header, b) < 0)
      wzfatal("Cannot write bam.\n");
	d->n_mapped++;

	return 0;
}

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: bsstrand [options] ref.fa in.bam [out.bam]\n");
  fprintf(stderr, "Input options:\n");
	fprintf(stderr, "     -g        region (optional, chrX:123-456 if missing, process the whole bam).\n");
	fprintf(stderr, "     -e        output abnormal read stats [optional].\n");
	fprintf(stderr, "     -a        output all reads stats [optional].\n");
	fprintf(stderr, "     -m        infer bsstrand information if missing [optional].\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

/* if nC2T > 1 and nG2A > 1: then fail and mark "?"
   if nC2T - nG2A > 2 and "-" => "+"
   if nG2A - nC2T > 2 and "+" => "-"
   output a summary of how many reads are inconsistent */
int main(int argc, char *argv[]) {
	int c;
	char *reg = 0;										/* region */
	bsstrand_conf_t conf = {.output_read = 0, .output_all_read = 0, .infer_bsstrand = 0};

  if (argc < 2) return usage();
  while ((c = getopt(argc, argv, "i:o:r:g:eamh")) >= 0) {
    switch (c) {
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
  char *reffn = optind < argc ? argv[optind++] : NULL;
  char *infn  = optind < argc ? argv[optind++] : NULL;
  char *outfn = optind < argc ? argv[optind++] : NULL;
  if (!reffn || !infn) {
    usage();
    wzfatal("Please provide reference and input bam.\n");
  }

	bsstrand_data_t d = {
		.n_corr = 0,
		.n_mapped = 0,
		.n_fail = 0,
		.n_unmapped = 0,
		.rs = init_refcache(reffn, 100, 100000),
		.conf = &conf,
	};
	bam_filter(infn, outfn, reg, &d, bsstrand_func);

  /*** output stats ***/
	fprintf(stderr, "Mapped reads: %d\n", d.n_mapped);
	fprintf(stderr, "Unmapped reads: %d\n", d.n_unmapped);
	fprintf(stderr, "Corrected reads: %d (%1.2f%%)\n", d.n_corr, (double)d.n_corr/(double)d.n_mapped*100.);
	fprintf(stderr, "Failed reads: %d (%1.2f%%)\n", d.n_fail, (double)d.n_fail/(double)d.n_mapped*100.);
	free_refcache(d.rs);
	return 0;
}
