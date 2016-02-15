#ifndef _WZ_REFSEQ_H_
#define _WZ_REFSEQ_H_

#include "faidx.h"

typedef struct {
	faidx_t *fai;
	char *chrm;
	uint32_t beg;
	uint32_t end;
	char *seq;
	uint32_t flank1;
	uint32_t flank2;
	int seqlen;
} refseq_t;

static inline refseq_t* init_refseq(char *ref_fn, uint32_t flank1, uint32_t flank2) {
	
	refseq_t *rs = calloc(1, sizeof(refseq_t));
  rs->fai = fai_load(ref_fn);
  if (!rs->fai) {
    fprintf(stderr, "[%s:%d] Cannot load reference %s\n", __func__, __LINE__, ref_fn);
    fflush(stderr);
    exit(1);
  }
	rs->flank1 = flank1;
	rs->flank2 = flank2;
	return rs;
}


/* beg and end are 1-based */
static inline void fetch_refseq(refseq_t *rs, char *chrm, uint32_t beg, uint32_t end) {

	if (rs->chrm != 0
			&& strcmp(chrm, rs->chrm) == 0
			&& rs->beg <= beg
			&& rs->end >= end) return;
	else {

		/* get sequence length */
		rs->seqlen = faidx_seqlen(rs->fai, chrm);
		if (rs->seqlen < 0) {
			fprintf(stderr, "[%s:%d] Error, cannot retrieve reference %s:%u-%u.\n",
							__func__, __LINE__, chrm, rs->beg, rs->end);
			exit(1);
		}

		/* beg and end */
		if (rs->flank1 > beg + 1) rs->beg = 1;
		else rs->beg = beg - rs->flank1;
		if (end + rs->flank2 > (unsigned) rs->seqlen) rs->end = rs->seqlen;
		else rs->end = end + rs->flank2;

		rs->chrm = realloc(rs->chrm, strlen(chrm)+1);
		strcpy(rs->chrm, chrm);
		if (rs->seq) free(rs->seq);
		int l;
		rs->seq = faidx_fetch_seq(rs->fai, rs->chrm, rs->beg-1, rs->end-1, &l);
		if ((unsigned) l != rs->end-rs->beg+1){
			fprintf(stderr, "[%s:%d] Error, cannot retrieve reference %s:%u-%u.\n",
							__func__, __LINE__, chrm, rs->beg, rs->end);
			exit(1);
		}
	}
}

static inline void free_refseq(refseq_t *rs) {

	if (rs->seq) free(rs->seq);
	if (rs->chrm) free(rs->chrm);
	fai_destroy(rs->fai);
  free(rs);

}

static inline int in_range_refseq(refseq_t *rs, uint32_t rpos) {
  if (rpos < rs->beg || rpos > rs->end) return 0;
  else return 1;
}

/* rpos is 1-based */
static inline char getbase_refseq(refseq_t *rs, uint32_t rpos) {
	if (rpos<rs->beg || rpos>rs->end) {
		fprintf(stderr, "[%s:%d] Error retrieving base %u outside range %s:%u-%u.\n",
						__func__, __LINE__, rpos, rs->chrm, rs->beg, rs->end);
		exit(1);
	}
	return rs->seq[rpos-rs->beg];
}

/* rpos is 1-based */
static inline char *subseq_refseq(refseq_t *rs, uint32_t rpos) {
	return rs->seq+rpos-rs->beg;
}

static inline void subseq_refseq2(refseq_t *rs, uint32_t rpos, char *seq, int len) {
	if (rpos < rs->beg || rpos+len-1 > rs->end) {
		fprintf(stderr, "[%s:%d] Error retrieving range %u-%u outside range %s:%u-%u.\n",
						__func__, __LINE__, rpos, rpos+len, rs->chrm, rs->beg, rs->end);
		exit(1);
	}
	int i;
	for (i=0; i<len; ++i)
		seq[i] = toupper(rs->seq[rpos-rs->beg+i]);
}

#endif /* _WZ_REFSEQ_H_ */
