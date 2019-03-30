/* Rectangularize an epiread output */

#include "wztsv.h"
#include "wvec.h"
#include "kstring.h"
#include "refcache.h"

int refcache_next_cg(refcache_t *rc, int pos) {
   for (;; pos++) {
      if (toupper(refcache_getbase_auto(rc, NULL, pos)) == 'C' &&
          toupper(refcache_getbase_auto(rc, NULL, pos+1)) == 'G')
         return pos;
   }
}

static int usage() {
   fprintf(stderr, "\n");
   fprintf(stderr, "Convert epiread into a rectanglular matrix.\n");
   fprintf(stderr, "Usage: biscuit rectangle [options] [ref.fa] [.epiread]\n");
   fprintf(stderr, "Options:\n\n");
   fprintf(stderr, "     -o        output file [stdout]\n");
   fprintf(stderr, "     -h        this help.\n");
   fprintf(stderr, "\n");
   return 1;
}

DEFINE_VECTOR(kstring_v, kstring_t)

int main_rectangle(int argc, char *argv[]) {

   char *out_fn = 0;
   
   if (argc<2) return usage();
   int c;
   while ((c = getopt(argc, argv, "o:"))>=0) {
      switch (c) {
      case 'o': out_fn = optarg; break;
      case 'h': return usage();
      default:
         fprintf(stderr, "Unrecognized command\n");
      }
   }

   char *ref_fn = argv[optind++];
   char *epiread_fn = argv[optind++];
   refcache_t *rc = init_refcache(ref_fn, 1000, 1000);
   tsv_t *tsv = tsv_open(epiread_fn);

   int region_beg = 0, region_width = -1;
   kstring_v *reads = init_kstring_v(1000);
   char *chrm = NULL;
   while (tsv_read(tsv)) {
      int read_beg = atoi(tsv->fields[4]);
      if (!region_beg) region_beg = read_beg;

      if (chrm == NULL) {
         chrm = calloc(strlen(tsv->fields[0]+1), sizeof(char));
         strcpy(chrm, tsv->fields[0]);
         refcache_set_chromosome(rc, chrm);
      } else if (strcmp(chrm, tsv->fields[0]) != 0) {
         fprintf(
            stderr,
            "[%s:%d] Error, rectangle cannot cross chromosomes.\n",
            __func__, __LINE__);
         exit(1);
      }

      // compute padding
      int p, pad = 0;
      for (p = region_beg; p < read_beg; ++pad)
         p = refcache_next_cg(rc, p) + 1;

      kstring_t *read = next_ref_kstring_v(reads);
      while(pad--) kputc('N', read);
      kputs(tsv->fields[5], read);

      if (region_width < 0 || (unsigned) region_width < read->l)
         region_width = read->l;
   }

   FILE *out;
   if (out_fn) out = fopen(out_fn, "w");
   else out = stdout;

   /* Pad the reads with length */
   unsigned i;
   for (i = 0; i < reads->size; ++i) {
      kstring_t *read = ref_kstring_v(reads, i);
      while (read->l < (unsigned) region_width) kputc('N', read);
      fputs(read->s, out);
      fputc('\n', out);
   }

   tsv_close(tsv);
   
   return 0;
}
