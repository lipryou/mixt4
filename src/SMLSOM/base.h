#ifndef INCLUDE_base_h_
#define INCLUDE_base_h_

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <float.h>


#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand() //ˆê—l—”

#define EPS 1e-4

int debug;
int silent;

typedef struct list {
  int node;
  struct list *next;
} list;
typedef list* link;

typedef enum {false, true} boolean;
typedef double (*distfunc)(double,double *);
typedef void (*updfunc)(double *, double *, double *, link, int, int, int);

#endif //INCLUDE_base_h_
