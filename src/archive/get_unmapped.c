#include <stdio.h>
#include <string.h>
#include "sam.h"


/* this utility gets unmapped reads from the end of bam file assuming existence of chrM */

void getunmapped(char *infn, char *reg) {
	samfile_t *in = samopen(infn, "rb", 0);
	samfile_t *out = samopen("-", "wh", in->header);
	bam_index_t *idx = bam_index_load(infn);
	int tid, beg, end;
	bam_parse_region(in->header, reg, &tid, &beg, &end);
	bam_iter_t iter = bam_iter_query(idx, tid, beg ,end);
	bam1_t *b = bam_init1();
	int ret;
	bam_iter_seek(in->x.bam, iter);
	int cnt_unmapped = 0;
	while ((ret = bam_read1(in->x.bam, b)) >= 0) {
		if (b->core.flag & BAM_FUNMAP) {
			samwrite(out, b);
			/* char *s = bam_format1_core(in->header, b, BAM_OFDEC); */
			/* puts(s); */
			/* free(s); */
			/* cnt_unmapped++; */
		}
	}
	bam_destroy1(b);
	bam_index_destroy(idx);
	samclose(in);
	samclose(out);
	fprintf(stderr, "[%s:%d] Found %d unmapped reads\n", __func__, __LINE__, cnt_unmapped);
	fflush(stderr);
}


static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: getunmapped [options] -i [in.bam] -g chrM\n");
  fprintf(stderr, "Input options:\n");
	fprintf(stderr, "     -i        input bam.\n");
	fprintf(stderr, "     -g        region (optional, jump to region before processing), typically a late chromosome.\n");
  fprintf(stderr, "     -h        this help.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {
	int c;
	char *reg = 0;										/* region */
	char *infn = 0;
  if (argc < 2) return usage();
  while ((c = getopt(argc, argv, "i:g:h")) >= 0) {
    switch (c) {
		case 'i': infn = optarg; break;
		case 'g': reg = optarg; break;
    case 'h': return usage();
    default:
      fprintf(stderr, "[%s:%d] Unrecognized command: %c.\n", __func__, __LINE__, c);
      exit(1);
      break;
    }
  }

	getunmapped(infn, reg);

	return 0;
}
