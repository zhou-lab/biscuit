#include <stdio.h>
#include "gsl/gsl_cdf.h"

int main(int argc, char *argv[])
{
  fprintf(stdout, "test chisq P: %1.3f\n", gsl_cdf_chisq_P(10.0,2.0));
  fprintf(stdout, "test chisq Q: %1.3f\n", gsl_cdf_chisq_Q(10.0,2));
  return 0;
}
