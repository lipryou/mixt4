#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

extern double **L;
extern double **diagL;

extern void init_params(double *Covs, int Kmax, int p);
extern void free_params(int M, int p);
extern double mvnorm(double *xi, double *mu, int index, int p);

int main() {
  int i;
  int Kmax = 2;
  int p = 2;
  double *Covs = (double *)malloc(Kmax*p*p*sizeof(double));
  Covs[0] = 4;   Covs[1] = 1;   Covs[2] = 1;   Covs[3] = 2;
  Covs[4] = 2;   Covs[5] = 0.5;   Covs[6] = 0.5;   Covs[7] = 1;

  double mu[2] = {2.0, 2.0};
  double xi[2] = {1.0, 1.5};

  init_params(Covs, Kmax, p);

  // normal
  for (i=0; i < Kmax; i++)
    printf("the %dth comp. loglik: %.3f\n", i, mvnorm(xi, mu, i, p));

  free(Covs);
  free_params(Kmax, p);

  return 0;
}
