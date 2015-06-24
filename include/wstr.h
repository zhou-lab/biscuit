
#ifndef _WZ_STR_H_
#define _WZ_STR_H_

#include <stdio.h>
#include <stdarg.h>
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


#endif /* _WZ_STR_H_ */
