#ifndef INCLUDE_update_h_
#define INCLUDE_update_h_

#include "base.h"

extern updfunc matchUpdate(char **dtype);
extern void Update_mom1(double *xi, double *moments, double *alpha, link neighs, int p, int K, int M);
extern void Update_mom2(double *xi, double *moments, double *alpha, link neighs, int p, int K, int M);
extern void Update_momK(double *xi, double *moments, double *alpha, link neighs, int p, int K, int M);
extern void Update_multinom(double *xi, double *theta, double *alpha, link neighs, int p, int K, int M);

#endif //INCLUDE_update_h_



