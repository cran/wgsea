W.combine <- function(W,n) {
  if(!is.list(W) || !is.numeric(n))
    stop("W must be a list of W values to combine, with n a vector of sample sizes.\n")
  if(length(W) != length(n))
    stop("need equal length W and n.\n")
  for(i in 1:length(W))
    W[[i]] <- W[[i]]/(n[[i]] + 1)
  if(length(W[[1]])>1)
    return(rowSums(do.call("cbind",W)))
  ##    W <- mapply(function(x,y) { x/y }, W, as.list(n))
  ##    W <- unlist(W)/n
  return(sum(unlist(W)))
}
Z.value <- function(W,Wstar,n.in,n.out) {
  if(is.list(W) && (!is.list(Wstar) || length(W)!=length(Wstar) || length(W)!=length(n.in) || length(W) !=length(n.out)))
    stop("W, Wstar must each be lists of length equal to the length of n.in and n.out.\n")
  if(is.list(W)) {
    W <- W.combine(W,n.in+n.out)
    Wstar <- W.combine(Wstar,n.in+n.out)
    mu.theor <- sum(mu.wilcox(n.in,n.out) / (n.in+n.out+1))
  } else {
    mu.theor <- mu.wilcox(n.in,n.out)
  }
  sd.est <- sd(Wstar)
  mu.est <- mean(Wstar)
  z <- (W-mu.theor)/sd.est
  Z.theoretical <- list(statistic=c(Z=z),
                        p.value=2*pnorm(abs(z),lower.tail=FALSE),
#                        parameter=c(mean.null=mu.theor),
                        method="Wilcoxon theoretical mean",
                        data.name=deparse(substitute(W)))
  z <- (W-mu.est)/sd.est
  Z.empirical <- list(statistic=c(Z=z),
                      p.value=2*pnorm(abs(z),lower.tail=FALSE),
#                      parameter=c(mean.null=mu.est),
                      method="Wilcoxon empirical mean",
                      data.name=deparse(substitute(W)))
  class(Z.theoretical) <- "htest"
  class(Z.empirical) <- "htest"
  return(list(Z.theoretical=Z.theoretical,
              Z.empirical=Z.empirical))
}
