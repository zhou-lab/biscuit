
#ifndef _WZ_STR_H_
#define _WZ_STR_H_

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include "kstring.h"

#define kstring_init(s) kstring_t (s); (s).s=0; (s).m=(s).l=0;

#define _WZ_STRLEN 100

static inline char*
wasprintf(const char *fmt, ...) {

  char *s = malloc(_WZ_STRLEN);	/* default length */

  va_list args;
  va_start(args, fmt);
  int l = vsnprintf(s, _WZ_STRLEN, fmt, args); /* if fail, l is the number characters that would've been written */
  if (l >= _WZ_STRLEN) {
    s = realloc(s, l+1);
    va_start(args, fmt);
    vsnprintf(s, l+1, fmt, args);
  }
  va_end(args);

  return s;
}

static inline int strcount_char(char *s, char c) {
  int i, n=0;
  for (i=0; s[i]; ++i)
    if (s[i] == c)
      ++n;
  return n;
}

static inline void ensure_number(char *s) {
  int i;
  for (i=0;s[i];++i) {
    if (!isdigit(s[i]) && s[i]!='.') {
      fprintf(stderr, "[%s:%d] Trying to convert nondigit string to number: %s\n", __func__, __LINE__, s);
      fflush(stderr);
      exit(1);
    }
  }
}

#endif /* _WZ_STR_H_ */
