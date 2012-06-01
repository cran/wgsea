\name{Z.value}
\alias{Z.value}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate a Z score from a Wilcoxon statistic and a set of random
  Wilcoxon statistics
}
\description{
  The mean of a Wilcoxon statistic is unaffected by correlation within
  the variable under test, but its variance is.  This function uses a
  set of Wilcoxon statistics generated from permuted data to estimate
  the variance empirically, and thus calculate a Z score.
}
\usage{
Z.value(W, Wstar, n.in, n.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
    Wilcoxon statistic for observed data.
}
  \item{Wstar}{
    A vector of Wilcoxon statistics for a set of permuted data.
}
  \item{n.in}{
    The number of items (SNPs) in the regions to be tested.
}
  \item{n.out}{
    The number of items (SNPs) in the control regions.
}
}
\value{
A list with two elements:
\item{Z.theoretical }{which uses the theoretical mean
  of the Wilcoxon distribution under the null generated from n.in, n.out above}
\item{Z.empirical }{which uses Wstar to calculate an empirical estimate
  of the mean of the Wilcoxon distribution under the null}
}
\note{
  The function can also deal with combining W statistics from multiple
  strata, as is typical in a meta analysis of GWAS data, using van Elteren's method.  Strata may be
  defined by different geography or different SNP chips.  
}
\author{
  Chris Wallace
}
\seealso{
  \code{\link{wilcoxon}}
}
\examples{
x <- exp(-rexp(1000)) # uniform
y <- exp(-rexp(1000,0.8)) # skewed towards 0
W <- wilcoxon(p=c(x,y),snps.in=1:1000)

p.perm <- matrix(sample(c(x,y),replace=TRUE,size=10000),ncol=5)
Wstar <- wilcoxon(p=p.perm,snps.in=1:1000)

Z.value(W=W, Wstar=Wstar, n.in=1000, n.out=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
