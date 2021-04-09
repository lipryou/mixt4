#ifndef INCLUDE_methodMapStr_h_
#define INCLUDE_methodMapStr_h_

#include "base.h"

extern int uDeath(double *data, double *moments, distfunc lld, int n, int p, int M, int K, int *cls, double *lliks, int *counts, double *consts, double mdl, char **dtype);
extern boolean uPrune(double *data, double *moments, distfunc lld, int n, int p, int M, int K, int *cls, double *consts, double beta);
extern int mvnDeath(double *data, double *mu1, double *mu2, int n, int p, int M, int *cls, double *lliks, int *counts, double mdl);
extern boolean mvnPrune(double *data, double *mu1, double *mu2, int n, int p, int M, int *cls, double beta);

#endif //INCLUDE_methodMapStr_h_
