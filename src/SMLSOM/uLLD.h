#ifndef INCLUDE_uLLD_h_
#define INCLUDE_uLLD_h_

#include "base.h"

extern double plldist(double,double *);
extern double glldist(double,double *);
extern double mltlldist(double f, double *pp);
extern double nblldist(double,double *);
extern double nlldist(double,double *);
extern double blldist(double f, double *pp);
extern distfunc match_dist(char **dtype);

#endif //INCLUDE_uLLD_h_
