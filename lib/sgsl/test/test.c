#include <stdio.h>
#include "gsl/gsl_cdf.h"

int main(int argc, char *argv[])
{
  fprintf(stdout, "testout: %1.3f\n", gsl_cdf_chisq_P(10.0,2.0));
  return 0;
}
