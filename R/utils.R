##初期値
LinearInit <- function(data, xdim=8, ydim=6) {
    pcm <- prcomp(data)
    pc <- pcm$rotation[, 1:2]
    sd <- pcm$sdev[1:2]
    mn <- apply(data, 2, mean)
    ans <- matrix(NA, xdim * ydim, dim(data)[2])
    colnames(ans) <- colnames(data)
    if (xdim >= ydim) {
        xtick <- sd[1] * pc[, 1]
        ytick <- sd[2] * pc[, 2]
    }
    else {
        xtick <- sd[2] * pc[, 2]
        ytick <- sd[1] * pc[, 1]
    }
    if (xdim == 1)
        xis <- rep(0, xdim)
    else xis <- seq(-2, 2, length = xdim)
    if (ydim == 1)
        yis <- rep(0, ydim)
    else yis <- seq(-2, 2, length = ydim)
    for (i in 1:(xdim * ydim)) {
        xi <- (i - 1)%%xdim + 1
        yi <- (i - 1)%/%xdim + 1
        ans[i, ] <- mn + xis[xi] * xtick + yis[yi] * ytick
    }
    ans
}

retK <- function(dtype) {
    switch(dtype,
           "pois" = 1,
           "geom" = 1,
           "nbin" = 2,
           "norm" = 2,
           "multinom" = 1,
           "binary" = 1
           )    
}

kmomcl <- function(data,class,k=1) {
    as.matrix(aggregate(data,by=list(class),function(x) mean(x^k))[,-1])
}

loglik <- function(data,class,thmat,dtype = c("pois","geom","nbin"),const) {
    dtype <- match.arg(dtype)
    logmat <- -llnorm(data,thmat,dtype)
    .const <- ifelse(dtype=="geom",0,const)
    sum(diag(as.matrix(aggregate(logmat,by=list(class),sum)[,-1]))) - .const
}

mdl <- function(data,class,thmat,dtype = c("pois","geom","nbin"),const) {
    dtype <- match.arg(dtype)
    K <- retK(dtype)
    N <- nrow(data)
    lik <- loglik(data,class,thmat,dtype,const)
    -lik + (length(thmat))/2 * log(N) + N*log(nrow(thmat)/K)
}
