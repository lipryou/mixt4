mixtures4 <- function(X,Kmax,Kmin=1,itmax=500,th=1e-4,silent=0) {
    MML <- function(llik,priors,k,n,p){
        -llik + dm*sum(log(priors)) + (dm + .5)*k*log(n)
    }
    ## only full covariance
    n <- nrow(X)
    p <- ncol(X)

    Mus <- X[sample(1:nrow(X),Kmax),]
    gCov <- cov(X)
    priors <- rep(1,Kmax)/Kmax

    Covs <- array(0,c(p,p,Kmax))
    for (k in 1:Kmax)
        Covs[,,k] <- diag(p) * max(diag(gCov/10))

    ## log-likelihood matrix k x n
    loglikmat <- log(priors) + t(mvnorm(X,Mus,Covs))
    ## likelihood matrix k x n
    weights <- exp(loglikmat)
    loglik <- sum(log(colSums(weights)))
    dm <- p*(p+3)/2
    mindl <- MML(loglik,priors,Kmax,n,p)

    dl <- c(mindl,rep(0,itmax-1))
    logliks <- c(loglik,rep(0,itmax-1))
    kappas <- c(Kmax,rep(0,itmax-1))
    ##browser()
    ret <- .C("mixtures",
              data = as.double(t(X)),
              weights = as.double(weights),
              Mus = as.double(t(Mus)),
              Covs = as.double(Covs),
              priors = as.double(priors),
              n = as.integer(n),
              p = as.integer(p),
              dmover2 = as.double(0.5*dm),
              K = as.integer(Kmax),
              Kmin = as.integer(Kmin),
              th = as.double(th),
              mindl = as.double(mindl),
              countf = integer(1),
              dl = as.double(dl),
              logliks = as.double(logliks),
              kappas = as.integer(kappas),
              trans1 = integer(itmax),
              trans2 = integer(itmax),
              blives = as.integer(rep(1, Kmax)),
              bMus = double(Kmax*p),
              bCovs = double(Kmax*p*p),
              bpriors = double(Kmax),
              itmax = as.integer(itmax),
              silent = as.integer(silent))

    countf <- ret$countf
    dl <- ret$dl[1:countf]
    logliks <- ret$logliks[1:countf]
    kappas <- ret$kappas[1:countf]
    trans1 <- which(ret$trans1 == 1)
    trans2 <- which(ret$trans2 == 1)
    blives <- ret$blives == 1
    Mus <- matrix(ret$bMus,Kmax,p,byrow=T)
    Covs <- array(ret$bCovs,c(p,p,Kmax))
    priors <- ret$bpriors[blives]
    Covs <- Covs[,,blives,drop=F]
    Mus <- Mus[blives,,drop=F]
    loglikmat <- log(priors) + t(mvnorm(X,Mus,Covs))
    weights <- exp(loglikmat)
    posteriors <- apply(weights,2,function(x) x/sum(x))
    if(is.vector(posteriors))
        posteriors <- matrix(posteriors,1,length(posteriors))
    clusters <- apply(posteriors,2,which.max)

    list(posteriors=posteriors,clusters=clusters,nc=sum(blives),Mus=Mus,Covs=Covs,priors=priors,dl=dl,logliks=logliks,kappas=kappas,trans1=trans1,trans2=trans2)
}

gmmEM <- function(X,Kmax,silent=0,itmax=500,th=1e-4) {
    gmmEM.BIC <- function(loglik,Kmax,n){
        -2*loglik + (Kmax-1 + Kmax*dm)*log(n)
    }
    gmmEM.ICL <- function(bic,posteriors) {
        C <- apply(posteriors,2,function(x) x==max(x))
        bic - 2*sum(C * ifelse(posteriors > 0,log(posteriors),0))
    }
    if (is.vector(X))
        X <- matrix(X,length(X),1)
    ## only full covariance
    n <- nrow(X)
    p <- ncol(X)

    Mus <- X[sample(1:nrow(X),Kmax),,drop=F]
    gCov <- cov(X)
    priors <- rep(1,Kmax)/Kmax

    Covs <- array(0,c(p,p,Kmax))
    for (k in 1:Kmax)
        Covs[,,k] <- diag(p) * max(diag(gCov/10))

    ## log-likelihood matrix k x n
    loglikmat <- log(priors) + t(mvnorm(X,Mus,Covs))
    ## likelihood matrix k x n
    weights <- exp(loglikmat)
    loglik <- sum(log(colSums(weights)))
    dm <- p*(p+3)/2

    logliks <- c(loglik,rep(0,itmax-1))
    ##browser()
    ret <- .C("EM_fullcovGMM",
              data = as.double(t(X)),
              weights = as.double(weights),
              Mus = as.double(t(Mus)),
              Covs = as.double(Covs),
              priors = as.double(priors),
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(Kmax),
              th = as.double(th),
              countf = integer(1),
              logliks = as.double(logliks),
              itmax = as.integer(itmax),
              silent = as.integer(silent))

    countf <- ret$countf
    logliks <- ret$logliks[1:countf]
    loglik <- tail(logliks,1)
    weights <- matrix(ret$weights,Kmax,n)
    Mus <- matrix(ret$Mus,Kmax,p,byrow=T)
    Covs <- array(ret$Covs,c(p,p,Kmax))
    priors <- ret$priors
    posteriors <- apply(weights,2,function(x) x/sum(x))
    if(is.vector(posteriors))
        posteriors <- matrix(posteriors,1,length(posteriors))
    bic <- gmmEM.BIC(loglik,Kmax,n)
    icl <- gmmEM.ICL(bic,posteriors)
    clusters <- apply(posteriors,2,which.max)

    list(posteriors=posteriors,clusters=clusters,Mus=Mus,Covs=Covs,priors=priors,logliks=logliks,loglik=loglik,bic=bic,icl=icl)
}

stepGmmEM <- function(X,K=1:9,seed=100,th=1e-4,itmax=500) {
    tmp.mat <- matrix(0,length(K),2)
    colnames(tmp.mat) <- c("BIC","ICL")
    for (k in K) {
        set.seed(seed+k)
        tmp <- gmmEM(X,k,silent=1,itmax=itmax,th=th)
        tmp.mat[k,"BIC"] <- tmp$bic
        tmp.mat[k,"ICL"] <- tmp$icl
    }
    tmp.mat
}

mvnorm <- function(X,Mus,Covs) {
    if (!is.array(Covs))
        stop("Covs should be an array")
    M <- nrow(Mus)
    p <- ncol(X)
    n <- nrow(X)
    ret <- matrix(0,n,M)
    for (m in 1:M) {
        tmp <- eigen(Covs[,,m])
        pos <- tmp$values > -10^-4
        tmp$values[pos] <- pmax(tmp$values[pos],10^-4)
        tmp$values[!pos] <- 0
        A <- apply(X,1,function(x) { xx <- t(tmp$vectors) %*% (x-Mus[m,]); (t(xx) %*% diag(tmp$values^-1) %*% xx) })
        ret[,m] <- -p/2*log(2*pi)-1/2*sum(log(tmp$values))-1/2 * A
    }
    ret
}

dyn.load("./mixt4.so")
