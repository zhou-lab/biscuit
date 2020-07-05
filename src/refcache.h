/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
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

#include <stdint.h>
#include "faidx.h"
#include "wzmisc.h"

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

static inline refcache_t* init_refcache(
   char *ref_fn, uint32_t flank1, uint32_t flank2) {
   
   refcache_t *rc = calloc(1, sizeof(refcache_t));
   rc->fai = fai_load(ref_fn);
   if (!rc->fai) {
      fprintf(stderr, "[%s:%d] Cannot load reference %s\n",
              __func__, __LINE__, ref_fn);
      fflush(stderr);
      exit(1);
   }
   rc->flank1 = flank1;
   rc->flank2 = flank2;
   return rc;
}

static inline void __refcache_fetch(refcache_t *rc) {
   if (rc->seq) free(rc->seq);
   int l;
   rc->seq = faidx_fetch_seq(
      rc->fai, rc->chrm, rc->beg-1, rc->end-1, &l);
   
   if ((unsigned) l != rc->end-rc->beg+1)
      wzfatal("[%s:%d] Error, cannot retrieve reference: %s:%u-%u.",
              __func__, __LINE__, rc->chrm, rc->beg, rc->end);
}

/* if [beg, end] is within [rc->beg, rc->end], do nothing
 * else fetch sequence from [beg - flank1, end + flank2]
 * beg and end are 1-based */
static inline void refcache_fetch(
   refcache_t *rc, char *chrm, uint32_t beg, uint32_t end) {

   if (rc->chrm != 0
       && strcmp(chrm, rc->chrm) == 0
       && rc->beg <= beg
       && rc->end >= end) return;

   /* get sequence length */
   rc->seqlen = faidx_seq_len(rc->fai, chrm);
   if (rc->seqlen < 0) {
      fprintf(
         stderr,
         "[%s:%d] Error, cannot retrieve reference %s:%u-%u.\n",
         __func__, __LINE__, chrm, rc->beg, rc->end);
      
      exit(1);
   }

   /* beg and end */
   if (rc->flank1 >= beg) rc->beg = 1;
   else rc->beg = beg - rc->flank1;
   if (end + rc->flank2 > (unsigned) rc->seqlen) rc->end = rc->seqlen;
   else rc->end = end + rc->flank2;

   if (chrm != rc->chrm) {
      rc->chrm = realloc(rc->chrm, strlen(chrm)+1);
      strcpy(rc->chrm, chrm);
   }
   
   __refcache_fetch(rc);
}

// set chromosome, set location to 1
#define refcache_set_chromosome(rc, chrm) refcache_fetch(rc,chrm,1,1)

// fetch the whole chromosome
static inline void refcache_fetch_chrm(refcache_t *rc, char *chrm) {
   rc->seqlen = faidx_seq_len(rc->fai, chrm);
   if  (rc->seqlen < 0)
      wzfatal("[%s:%d] Error, cannot retrieve chromosome: %s",
              __func__, __LINE__, chrm);
   
   if (rc->beg == 1 &&
       rc->end == (unsigned) rc->seqlen &&
       strcmp(chrm, rc->chrm) == 0) return;

   rc->beg = 1; rc->end = rc->seqlen;

   if (chrm != rc->chrm) {
      rc->chrm = realloc(rc->chrm, strlen(chrm)+1);
      strcpy(rc->chrm, chrm);
   }
   
   __refcache_fetch(rc);
}

static inline void free_refcache(refcache_t *rc) {

   if (rc->seq) free(rc->seq);
   if (rc->chrm) free(rc->chrm);
   fai_destroy(rc->fai);
   free(rc);

}


// [> see if the rpos is in range of rc <]
// static inline int in_range_refcache(refcache_t *rc, uint32_t rpos) {
//   if (rpos < rc->beg || rpos > rc->end) return 0;
//   else return 1;
// }


/************************************************/
/* The following does NOT auto-refetch sequence */
/************************************************/
/* rpos is 1-based */
static inline char refcache_getbase(refcache_t *rc, uint32_t pos) {
   
   if (pos<rc->beg || pos>rc->end)
      wzfatal(
         "[%s:%d] Error retrieving base %u outside range %s:%u-%u.\n",
         __func__, __LINE__, pos, rc->chrm, rc->beg, rc->end);
   
   return rc->seq[pos-rc->beg];
}

#define refcache_getbase_upcase(rc, pos) toupper(refcache_getbase(rc, pos))

/* rpos is 1-based */
static inline char *subseq_refcache(refcache_t *rc, uint32_t rpos) {
   return rc->seq+rpos-rc->beg;
}

/* get uppercased subsequence, checking range */
static inline void subseq_refcache2(
   refcache_t *rc, uint32_t rpos, char *seq, int len) {
   
   if (rpos < rc->beg || rpos+len-1 > rc->end) {
      fprintf(
         stderr,
         "[%s:%d] Error retrieving range %u-%u outside range %s:%u-%u.\n",
         __func__, __LINE__, rpos, rpos+len, rc->chrm, rc->beg, rc->end);
      exit(1);
   }
   int i;
   for (i=0; i<len; ++i)
      seq[i] = toupper(rc->seq[rpos-rc->beg+i]);
}


// if chrm is not given, then the function
// doesn't check chromosome, i.e.
// refcache_getbase_auto(rc, NULL, pos) assume
// taking the same chromosome
static inline char refcache_getbase_auto(
   refcache_t *rc, char *chrm, uint32_t pos) {

   if (chrm && strcmp(chrm, rc->chrm) != 0)
      refcache_fetch(rc, chrm, pos, pos);
   else if (pos < rc->beg || pos > rc->end)
      refcache_fetch(rc, rc->chrm, pos, pos);

   return refcache_getbase(rc, pos);
}

#endif /* _WZ_REFSEQ_H_ */
