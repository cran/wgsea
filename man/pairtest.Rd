\name{pairtest}
\alias{pairtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate p values for each SNP for case-control comparisons.
}
\description{
A wrapper for the \code{snpStats} function
\code{single.snp.tests}.   Generates p values for the association of each
SNP with case or control status.
}
\usage{
pairtest(case, control, n.perm = 0, pheno.perm = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{case}{
    \code{SnpMatrix} object holding genotypes of case subjects
}
  \item{control}{
    \code{SnpMatrix} object holding genotypes of control subjects
}
  \item{n.perm}{
    number of permutations of case control status required to generate
    permuted p value vectors.  The default, given by \code{n.perm=0}, is
    not to permute.
}
  \item{pheno.perm}{
    An alternative to specifying \code{n.perm} is to supply a matrix of
    alternative phenotypes, with each column relating to a different permutation.
}
}
\value{
  If \code{n.perm=0}, a vector of p values, one for each SNP (each
  column in the \code{case} and \code{control} objects.  If
  \code{n.perm>0}, a matrix of p values, each column representing the
  results of a different permutation.
 %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
  Chris Wallace
}
\examples{
data(for.exercise,package="snpStats")
case <- snps.10[subject.support$cc==1,]
control <- snps.10[subject.support$cc==0,]
summary(pairtest(case,control))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
