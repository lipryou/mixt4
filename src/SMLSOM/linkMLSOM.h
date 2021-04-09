#ifndef INCLUDE_linkMLSOM_h_
#define INCLUDE_linkMLSOM_h_

#include "base.h"

extern void ulinkMLSOM(double *data, double *moments, updfunc Update, distfunc lld, 
		       int *adjmatrix, double *alphas, double *radii, double *changes,
		       int n, int p, int M, int K, int rlen, int chgbyns, double *consts);
extern void mvn_linkMLSOM(double *data, double *mu1, double *mu2, int *adjmatrix,
			  double *alphas, double *radii, double *changes,
			  int n, int p, int M, int rlen, int chgbyns);

#endif //INCLUDE_linkMLSOM_h_
