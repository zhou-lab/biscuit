#include <stdio.h>
#include "wztsv.h"

int main(int argc, char *argv[])
{
  tsv_t *tsv = tsv_open(argv[1]);
  while (tsv_read(tsv)) {
    unsigned i;
    for (i=0; i < tsv->n; ++i) printf("\t>>|%s", tsv->fields[i]);
    printf("\n");
  }
  tsv_close(tsv);
  
  return 0;
}
