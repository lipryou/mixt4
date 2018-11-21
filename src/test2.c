#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

extern double **L;
extern double **diagL;

extern void init_params(double *Covs, int Kmax, int p);
extern void free_params(int M, int p);
extern int update_L(double *cov, int comp, int p);


int main() {
  int i, j, flg;
  int Kmax = 2;
  int p = 2;
  double *Covs = (double *)malloc(Kmax*p*p*sizeof(double));
  Covs[0] = 4;   Covs[1] = 1;   Covs[2] = 1;   Covs[3] = 2;
  Covs[4] = 2;   Covs[5] = 0.5;   Covs[6] = 0.5;   Covs[7] = 1;

  double mu[2] = {2.0, 2.0};
  double xi[2] = {1.0, 1.5};

  init_params(Covs, Kmax, p);

  // normal
  Covs[3] = 1;
  flg = update_L(Covs, 0, p);
  printf("successfuly? %d\n", flg);

  // abnormal
  Covs[3] = -1;
  flg = update_L(Covs, 0, p);
  printf("successfuly? %d\n", flg);

  free(Covs);
  free_params(Kmax, p);

  return 0;
}
