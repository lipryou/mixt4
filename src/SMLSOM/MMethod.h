#ifndef INCLUDE_MMethod_h_
#define INCLUDE_MMethod_h_

#include "base.h"

extern void ubestMatchs(double *data, double *moments, distfunc lld, int n, int p, int M, int K, 
			int *cls, double *lliks, int *counts, double *consts, link avails);
extern void ubatch(double *data, double *moments, distfunc lld, int n, int p, int M, int K, 
		   int *cls, double *lliks, int *counts, double *consts, link avails, char **dtype);
extern double uMDL(double *lliks, int n, int p, int K, int nnodes, char **dtype);
extern void mvn_bestMatchs(double *data, double *mu1, double *mu2, int n, int p, int M, int *cls, double *lliks, int *counts, link avails);  
extern void mvn_batch(double *data, double *mu1, double *mu2, int n, int p, int *cls, double *lliks, int *counts, link avails);
extern double mvnMDL(double *lliks, int n, int p, int nnodes);

#endif //INCLUDE_mvnMethod_h_
