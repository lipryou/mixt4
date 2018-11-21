/* C implemantation of Figueiredo, Mario A. T., and Anil K. Jain. "Unsupervised learning of finite mixture models." IEEE Transactions on pattern analysis and machine intelligence 24.3 (2002): 381-396.
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <float.h>

double ***L;
double **diagL;

double mvnorm(double *xi, double *mu, int index, int p);
int choldc(double **a, double *p, int n);
int update_L(double *cov, int comp, int p);
void init_params(double *Covs, int Kmax, int p);
void free_params(int M, int p);
double* estep_fullcov(double *weights, int n, int K, int Kmax);
double mstep_fullcov(double *data, int n, int p, int comp,
                     double *mu, double *cov, double *post,
                     int Kmax);
void lshift_vec(double *vec, int comp, int K);
void lshift_mat(double *mat, int comp, int n, int K, int Kmax);
void calc_weights(double *weights, double *priors,
                  double *Mus, int *lives,
                  double *data, int n, int p, int K, int Kmax);
double loglik_mvGMM(double *data, int n, int p,
                    double *Mus, double *weights, double *priors,
                    int *lives, int K, int Kmax);

// main
void mixtures(double *data, double *weights, double *Mus,
              double *Covs, double *priors, int *pn,
              int *pp, double *pdmover2, int *pK,
              int *pKmin, double *ptol, double *pmindl,
              int *pcountf, double *dl, double *logliks,
              int *kappas, int *trans1, int *trans2,
              int *lives, double *bMus, double *bCovs,
              double *bpriors, int *pitmax, int *pverbose)
{
  int j, k, l, m, comp, min_m, countf;
  int n = *pn, p = *pp, K = *pK, Kmax = *pK,
    Kmin = *pKmin, itmax = *pitmax, verbose = *pverbose;
  double *posts, *mu, *cov;
  double sum, tmp, dm, n_m;
  double mindl = *pmindl, loglik;
  double dmover2 = *pdmover2, tol = *ptol;
  int k_cont, killed;
  int *mapp, idx;

  //========debug=======//
  // print_debug(priors,weights,n,p,K);
  //========debug=======//

  mapp = (int *)malloc(Kmax * sizeof(int));
  for (k = 0; k < Kmax; k++) mapp[k] = k;

  // make variables for calculating likelihood
  init_params(Covs, Kmax, p);

  k_cont = 1;
  countf = 0;

  if (verbose != 0)
    printf("%03d : loglik = %.4f\n", countf, logliks[0]);
  while(k_cont) {
    do {
      idx = 0;
      while (idx < K) {
        comp = mapp[idx];

        // load mean and covariance of "comp"th component
        mu = &Mus[comp*p];
        cov = &Covs[comp*p*p];

        //E-step
        posts = estep_fullcov(weights, n, K, Kmax);
        //M-step
        n_m = mstep_fullcov(data, n, p, idx, mu, cov, posts, Kmax);

        // when priros are small, kill its components
        priors[idx] = (n_m - dmover2)/n;
        if (priors[idx] < 0)
          priors[idx] = 0.0;

        sum = 0.0;
        for (m = 0; m < K; m++)
          sum += priors[m];
        for (m = 0; m < K; m++)
          priors[m] /= sum;

        killed = 0;
        // elliminate the dead component
        if (priors[idx] == 0.0) {
          if (verbose != 0)
            printf("K = %d, comp = %d\n", K, mapp[idx]);

          trans1[countf] = 1;
          killed = 1;

          // left shift priors and posts
          lshift_vec(priors, idx, K);
          lshift_mat(posts, idx, n, K, Kmax);

          // weights update
          weights = posts;

          // off
          for (m = idx; m < K-1; m++) mapp[m] = mapp[m+1];
          lives[comp] = 0;

          // done remove
          K = K-1;
        }

        if (killed == 0) {
          calc_weights(weights, priors, Mus, lives, data, n, p, K, Kmax);
          idx++;
        }
      } //while(idx < K)
      countf++;

      // calc posterior for MML
      loglik = loglik_mvGMM(data, n, p, Mus, weights, priors, lives, K, Kmax);
      // save
      logliks[countf] = loglik;
      if (verbose != 0)
        printf("%03d : K=%2d, loglik = %.3f\n",countf,K,loglik);

      // calc Minimum Message Length
      sum = 0.0;
      for (m = 0; m < K; m++)
        sum += log(priors[m]);
      dl[countf] = -loglik + dmover2*sum + (dmover2 + 0.5)*K*log(n);
      kappas[countf] = K;

      // valid termination rule
      tmp = logliks[countf] - logliks[countf-1];
    } while(fabs(tmp/logliks[countf-1]) > tol);

    // save model if current is best
    if (dl[countf] < mindl) {
      mindl = dl[countf];
      for (m = 0; m < Kmax; m++) {
        bpriors[m] = priors[m];
        for (j = 0; j < p; j++) {
          bMus[m*p + j] = Mus[m*p + j];
          for (l = 0; l < p; l++)
            bCovs[m*p*p + j*p + l] = Covs[m*p*p + j*p + l];
        }
      }
    }

    // force to kill smallest component
    if (K > Kmin) {
      // find the smallest comp
      min_m = -1;
      dm = DBL_MAX;
      for (m = 0; m < K; m++) {
        if (dm > priors[m]) {
          min_m = m;
          dm = priors[m];
        }
      }
      if (min_m == -1)
        fprintf(stderr, "can't find the smallest component\n");

      // kill the comp
      comp = mapp[min_m];
      lives[comp] = 0;
      for (m = min_m; m < K-1; m++) mapp[m] = mapp[m+1];

      // flusleft
      lshift_vec(priors, min_m, K);
      lshift_mat(posts, min_m, n, K, Kmax);

      // done delete
      K = K-1;

      sum = 0.0;
      for (m = 0; m < K; m++)
        sum += priors[m];
      for (m = 0; m < K; m++)
        priors[m] /= sum;

      trans2[countf] = 1;

      //============debug===========//
      //print_debug(priors,weights,n,p,K);
      //============================//

      countf++;
      // calc posterior prob for MML
      loglik = loglik_mvGMM(data, n, p, Mus, weights, priors, lives, K, Kmax);
      logliks[countf] = loglik;
      // MML
      sum = 0.0;
      for (m = 0; m < K; m++)
        sum += log(priors[m]);
      dl[countf] = -loglik + dmover2*sum + (dmover2 + 0.5)*K*log(n);
      kappas[countf] = K;
    } else {
      k_cont = 0;
    }
    if (countf > (itmax-1))
      k_cont = 0;
  }

  *pcountf = countf;
  free_params(Kmax, p);
  free(mapp);
}


double mvnorm(double *xi, double *mu, int index, int p)
/*
  Calculate loglikelihood of multivariate normal distribution

  xi: ith p-dimensional data vector
  mu: the p-dimensional mean vector of the the 'index'th component
  index: specify the 'index'th component
  p: the number of dim
 */
{
  int i, j;
  double tmp;
  double lexp = 0.0;
  double det = 1.0;

  // calculate value of LLD
  // |A| = |LL^t| = |L|^2
  // |L| = prod(diag(L))
  // A was replaced in L^-1
  for (i = 0; i < p; i++) {
    tmp = 0.0;
    for (j = 0; j <= i; j++)
      tmp += L[index][i][j] * (xi[j] - mu[j]); //L^-1 (x-mu)
    lexp += tmp * tmp;
    det *= diagL[index][i];
  }

  return (-0.5*p*log(2*M_PI) - log(det) - 0.5*lexp);
}

int choldc(double **a, double *p, int n)
{
  int i, j, k;
  double sum;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      for (sum = a[i][j], k = i-1; k >= 0; k--) sum -= a[i][k] * a[j][k];
      if (i == j) {
        if (sum <= 0) {
          return -1;
        }
        p[i] = sqrt(sum);
      } else {
        a[j][i] = sum / p[i];
      }
    }
  }

  return 0;
}

int update_L(double *cov, int comp, int p)
/*
   cov : p-order covariance matrix
  comp : the number of component
     p : dimension
*/
{
  int i, j, k;
  int flg;
  double sum;

  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++)
      L[comp][i][j] = cov[i*p + j];
  }

  // cholesky decomposition of L, for calculating the likelihood
  // L = LL^t where L is a lower triangle matrix
  flg = choldc(L[comp], diagL[comp],p);
  if (flg == -1)
    fprintf(stderr, "covariance become singular\n");

  // calculate inverse of L
  // L was replaced in L
  for (i = 0; i < p; i++) {
    L[comp][i][i] = 1.0 / diagL[comp][i];
    for (j = i+1; j < p; j++) {
      sum = 0.0;
      for (k = i; k < j; k++)
        sum -= L[comp][j][k] * L[comp][k][i];
      L[comp][j][i] = sum / diagL[comp][j];
    }
  }

  return flg;
}

void init_params(double *Covs, int M, int p)
/*create triangular matrix*/
{
  int i, m;

  // working space for cholevsky decomposition
  L = (double ***)malloc(M * sizeof(double **));
  for (m = 0; m < M; m++) {
    L[m] = (double **)malloc(p * sizeof(double *));
    for (i = 0; i < p; i++)
      L[m][i] = (double *)malloc(p * sizeof(double));
  }

  // prepare diag
  diagL = (double **)malloc(M * sizeof(double *));
  for (m = 0; m < M; m++) {
    diagL[m] = (double *)malloc(p * sizeof(double));
    for (i = 0; i < p; i++)
      diagL[m][i] = 0.0;
  }

  // set variable L and diagL
  for (m = 0; m < M; m++)
    update_L(&Covs[m*p*p], m, p);

}

void free_params(int M, int p)
{
  int i,m;

  for (m = 0; m < M; m++) {
    free(diagL[m]);
    for (i = 0; i < p; i++)
      free(L[m][i]);


    free(L[m]);
  }
  free(L);free(diagL);
}

double* estep_fullcov(double *weights, int n, int K, int Kmax)
/* calculate posterior probabilities*/
{
  int i,k;
  double s;
  double *posts = weights;

  for (i = 0; i < n; i++) {
    s = 0.0;
    for (k = 0; k < K; k++)
      s += weights[i*Kmax + k];
    for (k = 0; k < K; k++)
      posts[i*Kmax + k] = weights[i*Kmax + k] / s;
  }

  return posts;
}

double mstep_fullcov(double *data, int n, int p,
                     int comp, double *mu, double *cov,
                     double *posts, int Kmax)
// return estimated sample size
{
  int i, j, l;
  double sum = 0.0;
  double tmp;

  // estimate mean and covariance
  //  preparation
  for (j = 0; j < p; j++) {
    mu[j] = 0.0;
    for (l = 0; l < p; l++)
      cov[j*p + l] = 0.0;
  }

  //  estimation for mean vector
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++)
      mu[j] += posts[i*Kmax + comp] * data[i*p + j];
  }
  //  estimate sample size
  for (i = 0; i < n; i++)
    sum += posts[i*Kmax + comp];
  //  normalize
  for (j = 0; j < p; j++)
    mu[j] /= sum;

  //  estimation covariance
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      tmp = posts[i*Kmax + comp] * (data[i*p + j] - mu[j]);
      for (l = 0; l < p; l++)
        cov[j*p + l] += tmp * data[i*p + l];
    }
  }
  //  normalize
  for (j = 0; j < p; j++) {
    for (l = 0; l < p; l++)
      cov[j*p + l] /= sum;
  }

  // update decomposed covariance
  update_L(cov, comp, p);

  return sum;
}

void lshift_vec(double *vec, int id, int K)
{
  int m;

  for (m = id; m < K-1; m++)
    vec[m] = vec[m+1];
  vec[m] = 0;
}

void lshift_mat(double *mat, int id, int n, int K, int Kmax)
{
  int i, m;

  for (i = 0; i < n; i++ ) {
    for (m = id; m < K-1; m++)
      mat[i*Kmax + m] = mat[i*Kmax + m+1];
    mat[i*Kmax + m] = 0;
  }
}

void calc_weights(double *weights, double *priors,
                  double *Mus, int *lives,
                  double *data, int n, int p, int K, int Kmax)
{
  int i, m;
  int k = -1;

  for (m = 0; m < K; m++) {
    for (k=k+1; lives[k]==0; k++);
    for (i = 0; i < n; i++)
      weights[i*Kmax+m] = priors[m] * exp(mvnorm(&data[i*p], &Mus[k*p], k, p));
  }
}

double loglik_mvGMM(double *data, int n, int p,
                    double *Mus, double *weights, double *priors,
                    int *lives, int K, int Kmax)
// return sum of multivariate normal loglikehood
{
  int i, m;
  int k = -1;
  double sum, llik, w;
  double loglik = 0.0;

  for (i = 0; i < n; i++) {
    sum = 0.0;
    for (m = 0; m < K; m++) {
      for (k=k+1; lives[k]==0; k++);
      llik = mvnorm(&data[i*p], &Mus[k*p], k, p);
      w =  priors[m] * exp(llik);
      weights[i*Kmax + m] = w;
      sum += w;
    }
    loglik += log(sum);
  }

  return loglik;
}
