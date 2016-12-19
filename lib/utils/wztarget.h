
typedef struct target_t {
  unsigned tid;
  char *name;
  uint32_t len;
} target_t;

DEFINE_VECTOR(target_v, target_t);

static inline destroy_target_v(target_v *targets) {
  unsigned i;
  for (i=0; i<targets->size; ++i) {
    target_t *t = ref_target_v(targets, i);
    if (t->name) free(t->name);
  }
  free_target_v(targets);
}

static inline unsigned locate_or_insert_target_v(target_v *targets, char *chr) {
  unsigned i;
  target_t *t;
  for (i=0; i<targets->size; ++i) {
    t = ref_target_v(targets, i);
    if (strcmp(t->name, chr)==0) return i;
  }
  t = next_ref_target_v(targets, i);
  t->name = strdup(chr);
  t->tid = targets->size-1;
  return t->tid;
}

static inline target_v *bamheader2targets(bam_header_t *header) {
  target_v *targets = init_target_v(50);
  target_t *t;
  int i;
  for (i=0; i<header->n_targets; ++i) {
    t = next_ref_target_v(targets);
    t->tid = i;
    t->name = in->header->target_name[i];
    t->len = in->header->target_len[i];
  }
  return targets;
}

static inline target_t *get_target(target_v *targets, char *chr) {
  unsigned i;
  target_t *t;
  for (i=0; i<targets->size; ++i) {
    t = ref_target_v(targets, i);
    if (strcmp(t->name, chr)==0) return t;
  }
  return NULL;
}
