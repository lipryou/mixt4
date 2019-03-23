from mixtures import mmlem

from sklearn.datasets import make_classification
from scipy.stats import multivariate_normal as mvnorm
import pandas as pd
import numpy as np


def logmvnorm(X, Mus, Covs):
    M = Mus.shape[0]
    n, p = X.shape
    logliks = np.zeros((n, M))
    for m in range(M):
        w, v = np.linalg.eig(Covs[m, :, :])
        w = np.where(w > 1e-4, w, 1e-4)
        Z = X - Mus[m, :]
        tmp = np.dot(Z, v).dot(np.diag(1/np.sqrt(w)))
        exps = (tmp*tmp).sum(axis=1)
        logliks[:, m] = -p/2*np.log(2*np.pi)-1/2*np.sum(np.log(w))-1/2*exps

    return logliks


def MML(llik, priors, n, K, dm):
    return -llik + dm*np.sum(np.log(priors)) + (dm + .5)*K*np.log(n)


def fit(data, Kmax=5, Kmin=1, itmax=500, th=1e-4, verpose=0):
    n, p = data.shape
    Mus = data.sample(Kmax).values
    gCov = data.cov().values
    priors = np.ones(Kmax)/Kmax
    X = data.values

    Covs = np.zeros(p*p*Kmax).reshape(Kmax, p, p)
    for k in range(Kmax):
        Covs[k, :, :] = np.identity(p) * max(np.diag(gCov)/10)

    # log-likelihood matrix n x K
    loglikmat = np.log(priors) + logmvnorm(X, Mus, Covs)

    # likelihoods n x K
    weights = np.exp(loglikmat)
    loglik = np.sum(np.log(weights.sum(axis=1)))
    dm = p*(p+3)/2

    mindl = MML(loglik, priors, n, Kmax, dm)

    # variables for the record
    dl = np.zeros(itmax)
    logliks = np.zeros(itmax)
    kappas = np.zeros(itmax)
    countf = np.int(0)
    trans1 = np.zeros(itmax, dtype=np.int)
    trans2 = np.zeros(itmax, dtype=np.int)
    lives = np.zeros(Kmax, dtype=np.int)
    bMus = np.zeros(Kmax*p)
    bCovs = np.zeros(Kmax*p*p)
    bpriors = np.zeros(Kmax)
    dl[0] = mindl
    logliks[0] = loglik
    kappas[0] = Kmax

    mmlem(X, weights, Mus, Covs, priors, n, p,
          dm/2, Kmax, Kmin, th, mindl, countf,
          dl, logliks, kappas, trans1, trans2,
          lives, bMus, bCovs, bpriors, itmax, 1)


if __name__ == "__main__":
    X = make_classification(n_samples=1000, n_features=2, n_informative=1,
                            n_redundant=0, n_repeated=0, n_classes=2,
                            n_clusters_per_class=1, weights=None,
                            flip_y=0.01, class_sep=1.0, hypercube=True,
                            shift=0, scale=1.0, shuffle=True,
                            random_state=1)
    data = pd.DataFrame(X[0], X[1])
    fit(data)
