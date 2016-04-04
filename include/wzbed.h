
#include "wztarget.h"

typedef struct {
  char *bed;
  gzFile FH;
  char *chrm;
  meth_obs1_t next;             /* next object */
  kstring_t nextchrom;          /* chromosome of next object */
} methbed_t;

static inline void bed_close(methbed_t *m) {
  gzclose(m->FH);
}

static inline void free_bed(methbed_t *m) {
  free(m->chrm);
  free(m->nextchrom.s);
  free(m);
}

typedef struct bed1_t {
  unsigned tid;
  int beg;
  int end;
  void *data;
} bed1_t;

DEFINE_VECTOR(bed1_v, bed1_t)

static inline int bed_parse1(char *line, target_v *targets, bed1_t *b, void (*dataparser)(bed1_t*, char**)) {

  char *tok; char **linerest;
  tok=strtok_r(line, "\t", linerest);
  b->tid = locate_or_insert_target_v(targets, tok);

  /* start */
  tok=strtok(NULL, "\t", linerest);
  ensure_number(tok);
  b->beg = atoi(tok);

  /* end */
  tok=strtok(NULL, "\t", linerest);
  ensure_number(tok);
  ob->pos = atoi(tok);

  if (dataparser)
    dataparser(b, linerest);
  else
    b->data = NULL;

  return 1;
}


static inline bed1_v *bed_read(char *bedfn) {

  target_v *targets = init_target_v(2);
  bed1_v *beds = init_bed_v(2);

  gzFile FH = gzopen(bedfn);
  kstring_t line;
  line.l = line.m = 0; line.s = 0;
  bed1_t *b;
  FILE *fh = open(argv[1],"r");
  while (1) {
    int c=gzgetc(FH);
    if (c=='\n' || c==EOF) {
      b = next_ref_bed_v(beds);
      bed_parse1(line.s, targets, b, NULL);
      line.l = 0;
      if (c==EOF) {
        break;
      }
    } else {
      kputc(c, &line);
    }
    free(line.s);
  }
  return beds;
}

bed1_v *target2bed(target_v *targets) {
  unsigned i;
  bed1_v *beds = init_bed1_v(2);
  for (i=0; i<targets->size; ++i) {
    target_t *t = ref_target_v(targets, i);
    bed1_t *b = next_ref_bed1_v(beds);
    b->tid = i;
    b->beg = 0;
    b->end = t->len;
    b->data = NULL;
  }
  return beds;
}

static inline void bamregion2bed(bed1_t *bed, target_v *targets, char *str) {

  char *s;
	int i, l, k, name_end;

	*ref_id = b->beg = b->end = -1;
	name_end = l = strlen(str);
	s = (char*)malloc(l+1);
	// remove space
	for (i = k = 0; i < l; ++i)
		if (!isspace(str[i])) s[k++] = str[i];
	s[k] = 0; l = k;
	// determine the sequence name
	for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
	if (i >= 0) name_end = i;
	if (name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for (i = name_end + 1; i < l; ++i) {
			if (s[i] == '-') ++n_hyphen;
			else if (!isdigit(s[i]) && s[i] != ',') break;
		}
		if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
		s[name_end] = 0;
    
		target_t *t = get_target(targets, s);
		if (!t) { // cannot find the sequence name
			t = get_target(targets, str); // try str as the name
			if (!t) {
        fprintf(stderr, "[%s:%d] fail to determine sequence name.\n", __func__, __LINE__);
        fflush(stderr);
				free(s); return -1;
			} else s[name_end] = ':', name_end = l;
		}
	} else t = get_target(targets, str);
  if (!t) {
    free(s); 
    return -1;
  }
  b->tid = t->tid;
	// parse the interval
	if (name_end < l) {
		for (i = k = name_end + 1; i < l; ++i)
			if (s[i] != ',') s[k++] = s[i];
		s[k] = 0;
		b->beg = atoi(s + name_end + 1);
		for (i = name_end + 1; i != k; ++i) if (s[i] == '-') break;
		b->end = i < k? atoi(s + i + 1) : 1<<29;
		if (b->beg > 0) b->beg--;
	} else b->beg = 0, b->end = 1<<29;
	free(s);
	return b->beg <= b->end? 0 : -1;
}

static inline bed1_v *bamregion2bedlist(target_v *targets, char *region) {
  bed1_v *beds = init_bed1_v(2);
  bed1_t *b = next_ref_bed1_v(beds);
  if (bam_region2bed(b, targets, region) < 0) {
    fprintf(stderr, "[%s:%d] failed to parse region\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  return beds;
}
