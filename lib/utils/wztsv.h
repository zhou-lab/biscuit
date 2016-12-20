DEFINE_VECTOR(charp_v, char*)

typedef struct tsv_t {
  charp_v *fields;
  FILE *fp;
  char *fn;
  int read_finished;
} tsv_t;

static inline tsv_t *tsv_open(char *fn) {
  tsv_t *tsv = calloc(1, sizeof(tsv_t));
  tsv->fn = fn;
  tsv->fp = gzopen(fn);
  tsv->fields = init_charp_v(10);
  return tsv;
}

static inline void tsv_close(tsv_t *tsv) {
  free_charp_v(tsv->fields);
  gzclose(tsv->fp);
  free(tsv);
}

static inline char *tsv_read(tsv_t *tsv) {
  if (read_finished) return NULL;
  if (tsv->n_fields
}


