\name{wilcoxon}
\alias{wilcoxon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Wilcoxon test statistic, with optional weights.
}
\description{
Calculate a Wilcoxon two group rank test statistic, with optional
propensity score weighting.
}
\usage{
wilcoxon(p, snps.in, weights = NULL, binsize = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
    a numeric vector of observed p values from a list of SNPs or a matrix, with each column
    representing a vector under a different permutation of the dataset.
    %%     ~~Describe \code{p} here~~
}
\item{snps.in}{
  a numeric vector indicating which members of p form the test group
  (their complement form the control group).
%%     ~~Describe \code{snps.in} here~~
}
  \item{weights}{
    optional propensity score weights.  These are binned according to
    binsize, and a weight calculated which is inversely proportional to the probability of
    sampling a member of the test group in that bin.    
    %%     ~~Describe \code{weights} here~~
}
  \item{binsize}{
    see weights, above.
}
}
\value{
  A numeric value or, if p is a matrix, a numeric vector.
}
\references{
  Propensity weights are described
  
  Rosenbaum, P. R. & Rubin, D. B. The central role of the propensity
score in observational studies for causal effects. Biometrika, 1983, 70,
41-55

Rosenbaum, P. R. Model-based direct adjustment. Journal of the American
Statistical Association, 1987, 82, 387-394
}
\author{
  Chris Wallace
}
\seealso{
  \code{\link{Z.value}}
}
\examples{
x <- exp(-rexp(1000)) # uniform
y <- exp(-rexp(1000,0.8)) # skewed towards 0
wilcoxon(p=c(x,y),snps.in=1:1000)

## note, should be equal to
wilcox.test(x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
