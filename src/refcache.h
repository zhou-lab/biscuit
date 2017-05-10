/**
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

#ifndef _WZ_REFSEQ_H_
#define _WZ_REFSEQ_H_

#include "faidx.h"

#define bscall(b, pos) seq_nt16_str[bam_seqi(bam_get_seq(b), pos)]

/**************
 * refcache_t *
 **************/

/* A local cache of reference sequence.
 * This avoid too many disk accesses. */
typedef struct {
  faidx_t *fai;
  char *chrm;
  uint32_t beg;
  uint32_t end;
  char *seq;
  uint32_t flank1;
  uint32_t flank2;
  int seqlen;
} refcache_t;

static inline refcache_t* init_refcache(char *ref_fn, uint32_t flank1, uint32_t flank2) {
  refcache_t *rs = calloc(1, sizeof(refcache_t));
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

/* if [beg, end] is within [rs->beg, rs->end], do nothing
 * else fetch sequence from [beg - flank1, end + flank2]
 * beg and end are 1-based */
static inline void fetch_refcache(refcache_t *rs, char *chrm, uint32_t beg, uint32_t end) {

  if (rs->chrm != 0
      && strcmp(chrm, rs->chrm) == 0
      && rs->beg <= beg
      && rs->end >= end) return;
  else {

    /* get sequence length */
    rs->seqlen = faidx_seq_len(rs->fai, chrm);
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

static inline void free_refcache(refcache_t *rs) {

  if (rs->seq) free(rs->seq);
  if (rs->chrm) free(rs->chrm);
  fai_destroy(rs->fai);
  free(rs);

}

/* see if the rpos is in range of rs */
static inline int in_range_refcache(refcache_t *rs, uint32_t rpos) {
  if (rpos < rs->beg || rpos > rs->end) return 0;
  else return 1;
}

/* rpos is 1-based */
static inline char getbase_refcache(refcache_t *rs, uint32_t rpos) {
  if (rpos<rs->beg || rpos>rs->end) {
    fprintf(stderr, "[%s:%d] Error retrieving base %u outside range %s:%u-%u.\n",
            __func__, __LINE__, rpos, rs->chrm, rs->beg, rs->end);
    exit(1);
  }
  return rs->seq[rpos-rs->beg];
}

/* rpos is 1-based */
static inline char *subseq_refcache(refcache_t *rs, uint32_t rpos) {
  return rs->seq+rpos-rs->beg;
}

/* get uppercased subsequence, checking range */
static inline void subseq_refcache2(refcache_t *rs, uint32_t rpos, char *seq, int len) {
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
