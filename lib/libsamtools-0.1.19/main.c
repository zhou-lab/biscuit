#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "sam.h"

#ifdef _TEST_BAMIO_
int main(int argc, char **argv) {

  samfile_t *in = samopen("001.bam", "rb", 0);
  bam_header_t *header=in->header;
  bam1_t *b=bam_init1();

  printf("number of targets: %d\n", header->n_targets);
  int j;
  for (j=0;j<header->n_targets;j++) {
    printf("target %d : %s (%d)\n", j, header->target_name[j], header->target_len[j]);
  }

  printf("l_text: %d, n_text: %d\n", header->l_text, header->n_text);
  printf("plain text: %s\n", header->text);

  bam_header_t new_header;
  new_header.n_targets=25;
  new_header.target_name = (char**) malloc(sizeof(char*) * 25);
  for (j=0; j<25; j++) {
    new_header.target_name[j] = (char *) malloc(100);
    strcpy(new_header.target_name[j], "something");
    /* new_header.target_len[j] = 100; */
  }

  /* new_header.dict=NULL; */
  /* new_header.hash=NULL; */
  /* new_header.rg2lib=NULL; */
  /* new_header.l_text=0; */
  /* new_header.n_text=0; */
  /* new_header.text=""; */
  /* printf("OK111\n"); */
  samfile_t *out = samopen("001c.bam", "wb", &new_header);

  while (samread(in, b) >= 0) {
    int i;
    bam1_t *d;
    const bam1_core_t *c = &d->core;
    c->tid=1;
    c->pos=22222;
    c->flag=94;
    c->qual=10;			/* mapping quality */
    c->l_qname=strlen(qname)+1;
    c->n_cigar=1;
    c->l_qseq=10;

    const char qname[]="abced";
    d->m_data=100;
    d->data_len = c->l_qname+c->n_cigar+c->l_qseq;
    d->data=(uint8_t*) malloc(d->m_data);
    strcpy(bam1_qname(d), qname);
    memcpy(bam1_cigar(d), 100M);
    memcpy(bam1_seq(d), "seq");
    memcpy(bam1_qual(d), "qua");
    memcpy(bam1_aux(d), "");
    samwrite(out, b);
    break;
    /* uint8_t *s = bam1_seq(b), *t = bam1_qual(b); */
    /* printf("tid: %s\n", in->header->target_name[c->tid]); */

    /* if ((c->flag & 0x0020) && (bam_cigar_oplen(bam1_cigar(b)[0]) == 76) && (c->isize < 500) && (c->isize > 0)) { */
    /*   /\* fprintf(stdout, "%s\t%d\t%d\t%d\t%" SCNu16 "\t%" SCNu16 "\n", bam1_qname(b), c->isize, c->pos, c->tid, 0x0010 & c->flag, 0x0020 \ */
    /* 	 & c->flag); *\/ */
    /*   /\* chromosome, position, insert size *\/ */
    /*   fprintf(stdout, "%s\t%d\t%d\n", (char*) header->target_name[c->tid], c->pos, c->isize); */
    /*   /\* for(i=0; i<c->n_cigar; i++) { *\/ */
    /*   /\*      printf("%d%c\t", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i])); *\/ */
    /*   /\* } *\/ */
    /*   /\* putc('\n', stdout); *\/ */
    /* } */

    /* break; */
    /* fwrite(bam1_qname(b), c->l_qname-1, sizeof(char),stdout); */
    /* fputc('\t',stdout); */
    /* for (i = 0; i < c->l_qseq; ++i) fputc(bam_nt16_rev_table[bam1_seqi(s, i)], stdout); */
    /* fputc('\t',stdout); */
    /* if (t[0] == 0xff) */
    /*   { */
    /*        fputs("*",stdout); */
    /*   } */
    /* else */
    /*   { */
    /*        for (i = 0; i < c->l_qseq; ++i) fputc(t[i] + 33,stdout); */
    /*   } */
    /* fputc('\n',stdout); */

  }

  samclose(in);
  samclose(out);
  bam_destroy1(b);

  /* bamFile bam_file = bam_open("001.bam", "r"); */
  /* bam_header_t *header; */

  /* if (bam_file==NULL) return -1; */
  /* if (b==NULL) return -1; */
  /* header=bam_header_read(bam_file); */
  /* while(bam_read1(bam_file, b) >= 0) { */
  /*   int i; */
  /*   const bam1_core_t *c = &b->core; */
  /*   uint8_t *s = bam1_seq(b), *t = bam1_qual(b); */
  /*   if ((c->flag & 0x0020) && (bam_cigar_oplen(bam1_cigar(b)[0]) == 76) && (c->isize < 500) && (c->isize > 0)) { */
  /*     /\* fprintf(stdout, "%s\t%d\t%d\t%d\t%" SCNu16 "\t%" SCNu16 "\n", bam1_qname(b), c->isize, c->pos, c->tid, 0x0010 & c->flag, 0x0020 \ */
  /* 	 & c->flag); *\/ */
  /*     /\* chromosome, position, insert size *\/ */
  /*     fprintf(stdout, "%s\t%d\t%d\n", (char*) header->target_name[c->tid], c->pos, c->isize); */
  /*     /\* for(i=0; i<c->n_cigar; i++) { *\/ */
  /*     /\*      printf("%d%c\t", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i])); *\/ */
  /*     /\* } *\/ */
  /*     /\* putc('\n', stdout); *\/ */
  /*   } */

  /*   /\* break; *\/ */
  /*   /\* fwrite(bam1_qname(b), c->l_qname-1, sizeof(char),stdout); *\/ */
  /*   /\* fputc('\t',stdout); *\/ */
  /*   /\* for (i = 0; i < c->l_qseq; ++i) fputc(bam_nt16_rev_table[bam1_seqi(s, i)], stdout); *\/ */
  /*   /\* fputc('\t',stdout); *\/ */
  /*   /\* if (t[0] == 0xff) *\/ */
  /*   /\*   { *\/ */
  /*   /\*        fputs("*",stdout); *\/ */
  /*   /\*   } *\/ */
  /*   /\* else *\/ */
  /*   /\*   { *\/ */
  /*   /\*        for (i = 0; i < c->l_qseq; ++i) fputc(t[i] + 33,stdout); *\/ */
  /*   /\*   } *\/ */
  /*   /\* fputc('\n',stdout); *\/ */
  /* } */
  /* bam_header_destroy(header); */
  /* bam_close(bam_file); */
  /* bam_destroy1(b); */

  return 0;
}

#endif
