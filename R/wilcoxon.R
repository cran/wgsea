mu.wilcox <-
function(n1,n2) {
  exp(log(n1)+log(n2)-log(2))
}
sd.wilcox <- function(n1,n2) {
  sqrt( n1*(n2) * (n1+n2+1)/12 )
}

calc.wilcoxon <-
function(p,snps.in,n0,w=NULL) {
  R <- rank(p)[snps.in]
  if(!is.null(w))
    R <- R * w
  return(sum(R) - n0 * (n0+1)/2)
}

wilcoxon <-
function(p,snps.in,weights=NULL,binsize=0.05) {
  n0 <- length(snps.in)
  if(!is.null(weights)) {
    n <- length(weights)
    bin <- cut(weights,seq(0,0.5,by=binsize))
    f <- table(bin)/n
    f0 <- table(bin[snps.in])/n0
    bin.in <- bin[snps.in]
    w <- f[bin.in]/f0[bin.in]
  } else {
    w <- NULL
  }
  
  if(is.matrix(p)) { # assume permutations in each column
    nc <- ncol(p)
    W <- numeric(nc)
    for(j in 1:nc) {
      R <- rank(p[,j])[snps.in]
      if(!is.null(weights))
        R <- R * w
      W[j] <- calc.wilcoxon(p[,j],snps.in,n0,w=w)
    }
  } else {
    W <- calc.wilcoxon(p,snps.in,n0,w=w)
  }
  return(W)
}
