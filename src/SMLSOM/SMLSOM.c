#include "base.h"
#include "link.h"
#include "searchGraph.h"
#include "methodMapStr.h"
#include "uLLD.h"
#include "update.h"
#include "mvnLLD.h"
#include "MMethod.h"
#include "linkMLSOM.h"

// Notice
// 1d-array containing n x p matrix as follow:
//   array = |1st-row|2nd-row|...|nth-row|
// e.g. array[i*n + j] : i,jth-element of the matrix

void uSMLSOM(double *data, double *moments, char **dtype, Sint *adjmatrix,
	     double *alphas, double *pbeta, double *radii, double *changes,
	     Sint *pn, Sint *pp, Sint *pM, Sint *pK, Sint *prlen, Sint *pchgbyns,
	     Sint *classes, double *logliks, double *consts, Sint *ptau, Sint *pdebug, Sint *psilent) 
{
  int n = *pn, p = *pp, M = *pM, K = *pK, rlen = *prlen, chgbyns = *pchgbyns;
  int t, tau = *ptau;
  int i,j;
  boolean ispr = false;
  double beta = *pbeta;

  distfunc lld = match_dist(dtype);
  updfunc Update = matchUpdate(dtype);

  debug = *pdebug;
  silent = *psilent;  

  // initialization of graph list
  init_list(M);
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      if (adjmatrix[i*M + j] == 1)
	createEdge(i,j);
    }
  }

  int *counts = (int *)calloc(M,sizeof(int));
  int nd, death;
  double mdl;

  RANDIN;
  
  t = 0;
  do {
    // update reference vectors
    ulinkMLSOM(data, moments, Update, lld, adjmatrix, alphas, radii, &changes[t*rlen*n/chgbyns], n, p, M, K, rlen, chgbyns, consts);
    ubestMatchs(data, moments, lld, n, p, M, K, classes, logliks, counts, consts, availNodes);
    ubatch(data, moments, lld, n, p, M, K, classes, logliks, counts, consts, availNodes, dtype);
    t++;    
    if (t == tau)
      break;
    // update map structure
    nd = getDim(availNodes);
    mdl = uMDL(logliks, n, p, K, nd, dtype);
    if (beta != -1)
      ispr = uPrune(data, moments, lld, n, p, M, K, classes, consts, beta);
    death = uDeath(data, moments, lld, n, p, M, K, classes, logliks, counts, consts, mdl, dtype);
  } while ((death != -1 || ispr) && nd > 1);

  //
  graph2adjmat(adjmatrix, M);
  
  free_list(M); free(counts);  

  RANDOUT;  
}

void mvn_SMLSOM(double *data, double *mu1, double *mu2, Sint *adjmatrix,
		double *alphas, double *pbeta, double *radii, double *changes,
		Sint *pn, Sint *pp, Sint *pM, Sint *prlen, Sint *pchgbyns,
		Sint *classes, double *logliks, Sint *ptau, Sint *pdebug, Sint *psilent)
// data is input  n x p matrix
// mu1,mu2 are first and second sample moment for multivariate gaussian. 
{
  int n = *pn, p = *pp, M = *pM, rlen = *prlen, chgbyns = *pchgbyns;
  int t, tau = *ptau;
  int m,i,j;
  boolean ispr = false;
  double beta = *pbeta;

  debug = *pdebug;
  silent = *psilent;  

  // initialization of paramaters for calculating likelihood
  init_params(M,p);

  // initialization of graph list
  init_list(M);
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      if (adjmatrix[i*M + j] == 1)
	createEdge(i,j);
    }
  }

  int *counts = (int *)calloc(M,sizeof(int));
  int nd, death;
  double mdl;

  RANDIN;

  // initialization of L. see mvnLLD.c
  for (m = 0; m < M; m++)
    flgs[m] = updateA(mu1,mu2,m,p);

  t = 0;
  do {
    t++;
    // update reference vectors
    mvn_linkMLSOM(data, mu1, mu2, adjmatrix, alphas, radii, changes, n, p, M, rlen, chgbyns); //近所のノードを更新 学習率なども radiiはRのradiusと同じ
    mvn_bestMatchs(data, mu1, mu2, n, p, M, classes, logliks, counts, availNodes); 
    mvn_batch(data, mu1, mu2, n, p, classes, logliks, counts, availNodes);
    if (t == tau)
      break;
    // update map structure
    nd = getDim(availNodes);
    mdl = mvnMDL(logliks, n, p, nd);
    if (beta != -1)
      ispr = mvnPrune(data,mu1,mu2,n,p,M,classes,beta);
    death = mvnDeath(data, mu1, mu2, n, p, M, classes, logliks, counts, mdl);
  } while ((death != -1 || ispr) && nd > 1);

  //
  free_params(M,p);
  //
  graph2adjmat(adjmatrix, M);
  
  free_list(M); free(counts);  

  RANDOUT;  
}
