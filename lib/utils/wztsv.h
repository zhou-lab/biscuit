/**
 * The MIT License (MIT)
 * Copyright (c) 2016-2017 Wanding.Zhou@vai.org
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
 * usage:
 * tsv_t *tsv = tsv_open("input.tsv");
 * while (tsv_read(tsv) != NULL) {
 *   for (i=0; i < tsv->n; ++i) printf("\t|%s", tsv->fields[i]);
 *   printf("\n");
 * }
 * tsv_close(tsv);
**/

#include <string.h>
#include <stdlib.h>
#include <zlib.h>

typedef struct tsv_t {
  char **fields;
  size_t n, m;                  /* number of fields, and allocated */
  gzFile fp;
  int finished;
  int line_max;
  char *line;
} tsv_t;

static inline tsv_t *tsv_open(const char *fn) {
  tsv_t *tsv = calloc(1, sizeof(tsv_t));
  tsv->fp = gzopen(fn, "r");
  tsv->fields = NULL;
  tsv->line_max = 10000;
  tsv->line = calloc(tsv->line_max, sizeof(char));
  tsv->finished = 0;
  return tsv;
}

/* tsv_close release the memory */
static inline void tsv_close(tsv_t *tsv) {
  unsigned i;
  for (i = 0; i < tsv->m; ++i)
    free(tsv->fields[i]);
  free(tsv->fields);
  free(tsv->line);
  gzclose(tsv->fp);
  free(tsv);
}

static inline void __tsv_clear_fields(tsv_t *tsv) {
  tsv->n = 0;
}

/* i is 0-based */
static inline void __tsv_set_field(tsv_t *tsv, int i, const char *s) {
  tsv->fields[i] = realloc(tsv->fields[i], strlen(s)+1);
  strcpy(tsv->fields[i], s);
}

static inline void __tsv_append_field(tsv_t *tsv, const char *s) {
  if (tsv->n == tsv->m) {
    tsv->fields = realloc(tsv->fields, (++tsv->m)*sizeof(char*));
    /* the newly incremented memory is always initialized */
    tsv->fields[tsv->m-1] = NULL;
  }
  __tsv_set_field(tsv, tsv->n++, s);
}

/* return NULL when done */
static inline int tsv_read(tsv_t *tsv) {
  if (tsv->finished)
    return 0;

  char *tok;
  char *line_null = tsv->line;
  *line_null = '\0';

  while (1) {
    int c=gzgetc(tsv->fp);
    if (c == '\n' || c == EOF) { /* process current line */

      __tsv_clear_fields(tsv);
      tok = strtok(tsv->line, "\t");
      while (tok != NULL) {
        __tsv_append_field(tsv, tok);
        tok = strtok(NULL, "\t");
      }
      if (c == EOF)
        tsv->finished = 1;
      break;
    } else {                    /* read line */
      if (line_null < tsv->line + tsv->line_max) {
        *line_null++ = c;
        *line_null = '\0';
      } else {
        tsv->line_max <<= 1;
        int line_len = line_null - tsv->line;
        tsv->line = realloc(tsv->line, tsv->line_max);
        line_null = tsv->line + line_len;
      }
    }
  }

  return 1;
}
