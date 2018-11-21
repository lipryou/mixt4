#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

extern double **L;
extern double **diagL;

extern void init_params(double *Covs, int Kmax, int p);
extern void free_params(int M, int p);
extern void calc_weights(double *weights, double *priors,
                         double *Mus, int *lives,
                         double *data, int n, int p, int K, int Kmax);
double* estep_fullcov(double *weights, int n, int K, int Kmax);

int main() {
  int i, j;
  int Kmax = 2;
  int K = 2;
  int p = 2;
  int n = 3;
  int lives[2] = {1, 1};
  double *posts;

  double *Covs = (double *)malloc(Kmax*p*p*sizeof(double));
  double *Mus = (double *)malloc(Kmax*p*sizeof(double));
  double *data = (double *)malloc(n*p*sizeof(double));
  double *weights = (double *)malloc(n*Kmax*sizeof(double));
  double *priors = (double *)malloc(Kmax*sizeof(double));


  Covs[0] = 4;   Covs[1] = 1;   Covs[2] = 1;   Covs[3] = 2;
  Covs[4] = 2;   Covs[5] = 0.5;   Covs[6] = 0.5;   Covs[7] = 1;

  Mus[0] = 2; Mus[1] = 2; Mus[2] = 3; Mus[3] = 4;

  priors[0] = 0.4; priors[1] = 0.6;

  data[0] = 2; data[1] = 3; data[2] = 1; data[3] = 1;
  data[4] = 3; data[5] = 2;

  init_params(Covs, Kmax, p);

  calc_weights(weights, priors, Mus, lives, data, n, p, K, Kmax);
  for (i = 0; i < n; i++) {
    for (j = 0; j < K; j++)
      printf("%.3f ", weights[i*K + j]);
    printf("\n");
  }

  posts = estep_fullcov(weights, n, K, Kmax);
  for (i = 0; i < n; i++) {
    for (j = 0; j < K; j++)
      printf("%.3f ", posts[i*K + j]);
    printf("\n");
  }


  free(Covs); free(Mus); free(data);
  free(weights); free(priors);
  free_params(Kmax, p);

  return 0;
}
