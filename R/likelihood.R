lld.mvnorm <- function(data,mu1,Sigs) {
    M <- nrow(mu1)
    p <- ncol(data)
    n <- nrow(data)
    ret <- matrix(0,n,M)
    for (m in 1:M) {
        tmp <- eigen(Sigs[[m]])
        pos <- tmp$values > -10^-4
        tmp$values[pos] <- pmax(tmp$values[pos],10^-4)
        tmp$values[!pos] <- 0
        A <- apply(data,1,function(x) { xx <- t(tmp$vectors) %*% (x-mu1[m,]); (t(xx) %*% diag(tmp$values^-1) %*% xx) })
        ret[,m] <- -p/2*log(2*pi)-1/2*sum(log(tmp$values))-1/2 * A
    }
    ret
}

llnorm <- function(fmat,thmat,dtype=c("pois","geom","nbin","norm","multinom")) {
    ##return log-likelihood dissimilarity between sample and node.
    a2mat <- function(x) {
        if (!is.matrix(x)) {
            if(is.vector(x))
                x <- matrix(x,nrow=1)        
        } else {
            x <- as.matrix(x)
        }
        x
    }
    if(missing(dtype))
        stop("select dtype")
    dtype <- match.arg(dtype)

    K <- retK(dtype)
    
    if (!is.numeric(fmat) | !is.numeric(thmat))
        stop("argument is not numeric")

    fmat <- a2mat(fmat)
    thmat <- a2mat(thmat)
    
    fnr <- nrow(fmat)
    thnr <- nrow(thmat)
    M <- thnr / K
    
    if (ncol(fmat) != ncol(thmat))
        stop("argument error")

    if (dtype != "norm") {
        if (any(fmat < 0) | any(thmat < -1e-32))
            stop("argument is irregular")        
        thmat <- pmax(thmat,1e-100)
        if (dtype == "pois")
            return (t(log(thmat) %*% t(fmat) - apply(thmat,1,sum)) - rowSums(lfactorial(fmat)))
        else if (dtype == "geom")
            return (t(-apply(log1p(thmat),1,sum) + (log(thmat) - log1p(thmat)) %*% t(fmat)))
        else if (dtype == "multinom")
            return (lgamma(rowSums(fmat)+1) - rowSums(lgamma(fmat+1)) + fmat %*% t(log(thmat)))
        else if (dtype == "nbin") {
            mu1 <- thmat[(1:M)*K-1,,drop=F]
            mu2 <- thmat[(1:M)*K,,drop=F]
            p <- mu1/(mu2-mu1^2)
            p <- pmin(p,1.0-10^-4)
            p <- pmax(p,10^-4)
            r <- (mu1*p)/(1-p)
            lmat <- apply(r,1,function(x) colSums(lgamma(t(fmat)+x)))
            if (is.vector(lmat))
                dim(lmat) <- c(1,M)
            return (t(t(lmat) - apply(lgamma(r),1,sum) + apply(r*log(p),1,sum) + log(1-p) %*% t(fmat)) - rowSums(lgamma(fmat+1)))
        }
    }
    else if (dtype == "norm") {
        mu1 <- thmat[(1:M)*K-1,,drop=F]
        mu2 <- thmat[(1:M)*K,,drop=F]
        sig <- mu2 - mu1^2
        pos <- sig > -10^-4
        sig[pos] <- pmax(sig[pos],10^-4)        
        ret <- matrix(0,fnr,M)
        for (m in 1:M)
            ret[,m] <- ncol(fmat)/2 * log(2*pi) -sum(1/2 * log(sig[m,])) - apply( (t(fmat)-mu1[m,])^2 / (2*sig[m,]),2,sum)
        return (ret)
    }
}

loglikmat <- function(X,sobj) {
    sclass <- class(sobj)
    if (sclass == "usmlsom") {
        K <- sobj$K
        nd <- sobj$nd    
        M <- sobj$M
        np <- ncol(sobj$codes)
        if (K == 1) {
            codes <- with(sobj,codes[lives,])
        } else if (K == 2) {
            codes <- matrix(0,2*nd,np)
            mu1 <- with(sobj,codes[2*(1:M)-1,][lives,])
            rownames(mu1) <- which(lives)
            mu2 <- with(sobj,codes[2*(1:M),][lives,])
            rownames(mu2) <- which(lives)
            codes[2*(1:nd)-1,] <- mu1
            codes[2*(1:nd),] <- mu2
        }                    
        else if (sclass == "smlsom.mvnorm") {
            mu1 <- with(sobj,mu1[lives,])
            mu2 <- with(sobj,mu2[lives])
            nd <- sobj$nd
            Sigs <- list()
            for (m in 1:nd)
                Sigs[[m]] <- mu2[[m]] - mu1[m,] %*% t(mu1[m,])
            return (lld.mvnorm(X,mu1,Sigs))
        } else {
            stop("error")
        }        
    }
    llnorm(X,codes,sobj$dtype)
}
