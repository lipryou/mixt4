#ifndef INCLUDE_mvnLLD_h_
#define INCLUDE_mvnLLD_h_

#include "base.h"

double ***A;
double **diagL;
int *flgs;

extern void init_params(int M, int p);
extern void free_params(int M, int p);
extern double lld_mvnorm(double *xi, double *mu1, int m, int p);
/* calculate LLD based on the negative loglikelihood of multivariate normal distribution
   LLD = 1/2 * log(det(Sigma)) + 1/2 * (x-mu)^t Sigma^-1 (x-mu)
   note MVNORM :
 f(x;mu,Sigma) = 1/((2*pi)^(p/2)) * 1/det(Sigma)^(1/2) * exp(-1/2 * (x-mu)^t * Sigma^-1 * (x-mu))
*/
// require A store L^-1 and diagL store diagonals of L
extern int choldc(double **a, double *p, int n);
// cholesky decomposition
// A = LL^t where L is a lower triangular matrix
// if A is not postive definite then return -1, otherwise return 0
// result values are stored in input a and p

extern int updateA(double *mu1, double *mu2, int m, int p);

#endif //INCLUDE_mvnLLD_h_
