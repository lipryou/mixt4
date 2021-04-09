library("kohonen")
library("gtools")
source("~/Document/Rfuncs/utils.R")
source("~/Document/Rfuncs/likelihood.R")

smlsom.mvn <- function(data, xdim=8, ydim=6, topo="h", grid = somgrid(xdim,ydim,topo),
                       rlen = 1, alpha = c(0.05, 0.01),beta=15,radius = NA, init1, init2,
                       init.type=c("linear","random"),chgbyns=nrow(data), tau=xdim*ydim, debug=F, silent=0)
{
    if (!is.numeric(data))
        stop("Argument data should be numeric")
    if(!(is.matrix(data) | is.table(data)))
        stop("Argument data should be matrix or table")
    n <- nrow(data)
    p <- ncol(data)    
    M <- xdim*ydim
    if (missing(init1)) {
        init.type <- match.arg(init.type)
        if (init.type == "random") {
            cl <- sample(1:M,n,replace=T)
            init1 <- as.matrix(aggregate(data,by=list(cl),mean)[,-1])
            mu1 <- init1
        } else if (init.type == "linear") {
            init1 <- LinearInit(data,xdim,ydim)
            mu1 <- init1
        }
    } else {
        if (is.matrix(init1))
            stop("Argument init1 should be matrix")
        if (nrow(init1) != M)
            stop("The number of rows of init1 is not equal to M")
        if (ncol(init1) != p)
            stop("The number of columns of init1 is not equal to p")        
        mu1 <- init1
    }
    if (missing(init2)) {
        mu2 <- list()
        for (m in 1:M)
            mu2[[m]] <- diag(1,p) + mu1[m,] %*% t(mu1[m,])
    } else {
        if (is.matrix(init2))
            stop("Argument init2 should be matrix")
        if (nrow(init2) != M)
            stop("The number of rows of init2 is not equal to M")
        if (ncol(init2) != p)
            stop("The number of columns of init2 is not equal to p")
        mu2 <- list()
        for (m in 1:M)
            mu2[[m]] <- diag(init2[m,])
    }

    nhbrdist <- as.matrix(dist(grid$pts))
    adjmatrix <- round(nhbrdist,1) == 1
    
    if (length(radius) == 1) {
        if (M <= 4)
            radius <- c(2,-2)
        else {
            if (is.na(radius)) {
                radius <- quantile(nhbrdist,.67) * c(1, -1)
            } else {            
                radius <- sort(radius * c(1,-1),decreasing=TRUE)
            }
        }
    }
    
    changes <- rep(0, ceiling(n/chgbyns * rlen))
    classes <- rep(0,n)
    logliks <- rep(0,n)
    res <- .C("mvn_SMLSOM",
              data = as.double(t(data)),
              mu1 = as.double(t(mu1)),
              mu2 = as.double(sapply(mu2,as.vector)),
              adjmatrix = as.integer(adjmatrix),
              alphas = as.double(alpha),
              beta = as.double(beta),
              radii = as.double(radius),
              changes = as.double(changes),
              n = as.integer(n),
              p = as.integer(p),
              M = as.integer(M),
              rlen = as.integer(rlen),
              chgbyns = as.integer(chgbyns),
              classes = as.integer(classes),
              logliks = as.double(logliks),
              tau = as.integer(tau),
              debug = as.integer(debug),
              silent = as.integer(silent))
    
    changes <- matrix(res$changes, ncol=1)
    adjmatrix = matrix(res$adjmatrix,nrow=M,byrow=T)
    mu1 <- res$mu1
    mu1 <- matrix(mu1,M,p,byrow=T)
    rownames(mu1) <- 1:M
    mu2 <- split(res$mu2,rep(1:M,rep(p*p,M)))
    mu2 <- lapply(mu2,function(x) matrix(x,p,p))
    colnames(mu1) <- colnames(init1)
    unit.classif <- factor(res$classes+1,levels=1:M)
    lives <- adjmatrix[,1] != -1
    names(lives) <- 1:M
    logliks = res$logliks
    nd <- sum(lives)
    df <- nd*p*(p+3)/2
    mdl <- -sum(logliks) + df/2 * log(n) + n*log(nd)    

    structure(list(grid = grid, mu1 = mu1, mu2 = mu2, changes = changes, adjmatrix = adjmatrix, lives = lives, alpha = alpha, radius = radius, unit.classif=unit.classif, logliks = logliks, mdl=mdl,nd=nd,M=M),
              class = "smlsom.mvnorm")
}


kmeans.lld <- function(X,sobj) {
    logcdmat <- loglikmat(X,sobj)
    convs <- sum(apply(logcdmat,1,min))
    cat(sprintf("%.2f",convs))
    while(T) {
        cl <- apply(logcdmat,1,which.min)        
        mu1 <- as.matrix(aggregate(X,by=list(cl),mean)[,-1])
        Sigs <- by(X,cl,var)
        logcdmat <- lld.mvnorm(X,mu1,Sigs)    
        convs <- c(convs,-sum(apply(logcdmat,1,min)))
        cat(sprintf("%.2f",tail(convs,1)))
        tmp <- tail(convs,2)        
        if (abs((tmp[2]-tmp[1])/tmp[1]) < 10^-4)
            break
    }
    cl <- apply(logcdmat,1,which.min)
    structure(list(cl=cl,convs=convs,mu1=mu1,Sigs=Sigs))
}

usmlsom <- function(data, xdim=8, ydim=6, topo="h", grid = somgrid(xdim,ydim,topo),
                    rlen = 1, alpha = c(0.05, 0.01), beta = 15, dtype = c("pois","geom","nbin","norm","multinom", "binary"),
                    radius = NA, init1, init.type=c("linear","random"),
                    chgbyns=nrow(data), consts, tau=xdim*ydim, debug=F, silent=0)
{
    if (missing(dtype))
        stop("select dtype")
    dtype <- match.arg(dtype)
    K <- retK(dtype)

    ##const for calculate log-likelihood
    if (missing(consts)) {
        consts <- switch(dtype,
                         "pois" = apply(lfactorial(data),1,sum),
                         "geom" = rep(0,nrow(data)),
                         "nbin" = apply(lgamma(data+1),1,sum),
                         "norm" = rep(ncol(data)/2 * log(2*pi),nrow(data)),
                         "multinom" = rowSums(lgamma(data+1)) - lgamma(rowSums(data)+1),
                         "binary" = rep(0,nrow(data))
                         )
    }
    
    if (!is.numeric(data))
        stop("Argument data should be numeric")
    if(!(is.matrix(data) | is.table(data)))
        stop("Argument data should be matrix or table")
    
    n <- nrow(data)
    p <- ncol(data)    
    M <- xdim*ydim

    codes <- matrix(0,M*K,p)   
    if (missing(init1)) {
        if (dtype == "multinom") {
            init1 <- rdirichlet(M,rep(1,p))
        } else {
            init.type <- match.arg(init.type)
            if (init.type == "random") {
                cl <- sample(1:M,n,replace=T)
                init1 <- as.matrix(aggregate(data,by=list(cl),mean)[,-1])
                if (K >= 2) {
                    for (k in 2:K)
                        codes[0:(M-1)*K+k,] <- as.matrix(aggregate(data,by=list(cl),function(x) mean(x^k))[,-1])
                }
            } else if (init.type == "linear") {
                init1 <- LinearInit(data,xdim,ydim)
                initk <- sapply(1:K,function(k) mean(data^k))
                if (K >= 2) {
                    for (k in 2:K)
                        codes[0:(M-1)*K+k,] <- initk[k]
                }
            }
        }
    }
    codes[0:(M-1)*K+1,] <- init1

    colnames(codes) <- colnames(init1)

    nhbrdist <- as.matrix(dist(grid$pts))
    adjmatrix <- round(nhbrdist,1) == 1

    if (length(radius) == 1) {
        if (M <= 4)
            radius <- c(2,-2)
        else {
            if (is.na(radius)) {
                radius <- quantile(nhbrdist,.67) * c(1, -1)
            } else {            
                radius <- sort(radius * c(1,-1),decreasing=TRUE)
            }
        }
    }

    if (is.vector(alpha))
        alpha = matrix(alpha,2,K)
    if (is.matrix(alpha)) {
        if (nrow(alpha) != 2 || ncol(alpha) != K)
            stop("Arugement error : alpha")
    } else {
        stop("Argument alpha should be vector or matrix")
    }
           
    changes <- rep(0, tau*ceiling(n/chgbyns * rlen))
    classes <- rep(0,n)
    logliks <- rep(0,n)    
    res <- .C("uSMLSOM",
              data = as.double(t(data)),
              codes = as.double(codes),
              dtype = as.character(dtype),
              adjmatrix = as.integer(adjmatrix),              
              alphas = as.double(alpha),
              beta = as.double(beta),
              radii = as.double(radius),
              changes = as.double(changes),
              n = as.integer(n),
              p = as.integer(p),
              M = as.integer(M),
              K = as.integer(K),
              rlen = as.integer(rlen),
              chgbyns = as.integer(chgbyns),              
              classes = as.integer(classes),
              logliks = as.double(logliks),
              consts = as.double(consts),   
              tau = as.integer(tau),
              debug = as.integer(debug),
              silent = as.integer(silent))

    adjmatrix = matrix(res$adjmatrix,nrow=M,byrow=T)    
    codes <- res$codes
    dim(codes) <- c(M*K,p)
    colnames(codes) <- colnames(init1)
    unit.classif <- factor(res$classes+1,levels=1:M)    
    lives <- adjmatrix[,1] != -1
    names(lives) <- 1:M
    logliks = res$logliks
    nd <- sum(lives)
    changes <- matrix(res$changes[1:(ceiling(rlen*n/chgbyns)*(M-nd+1))], ncol=1)
    if (dtype == "multinom")
        df <- nd*(p-1)
    else
        df <- nd*K*p
    mdl <- -sum(logliks) + df/2 * log(n) + n*log(nd)    
    
    structure(list(grid = grid, codes = codes, changes = changes, adjmatrix = adjmatrix, lives = lives, alpha = alpha, radius = radius, unit.classif=unit.classif, logliks = logliks, mdl=mdl, nd=nd, K=K, M=M, dtype=dtype),
              class = "usmlsom")
}

##batchSOM <- function (data, xdim, ydim, dtype=c("pois","geom","nbin","norm","multinom"))
##{
##    grid = somgrid(xdim,ydim,"hex")
##    data <- as.matrix(data)
##    n <- nrow(data)
##    p <- ncol(data)
##    M <- nrow(grid$pts)    
##
##    dtype <- match.arg(dtype)
##    K <- retK(dtype)
##
##    thmat <- matrix(0,K*M,p)
##    thmat[K*(0:(M-1))+1,] <- LinearInit(data,xdim,ydim)        
##    if (K > 1)
##        thmat[K*(0:(M-1))+2,] <- matrix(colMeans(data^2),nrow=M,byrow=T)
##    nhbrdist <- as.matrix(dist(grid$pts))
##    for (r in radii) {
##        cl <- apply(loglikmat(X
##        A <- (nhbrdist <= r)[, cl]
##        ind <- rowSums(A) > 0
##        init[ind, ] <- A[ind, ] %*% data/rowSums(A)[ind]
##    }
##    structure(list(grid = grid, codes = init), class = "SOM")
##}

dyn.load("~/Document/Rfuncs/src/SMLSOM/SMLSOM.dll")
