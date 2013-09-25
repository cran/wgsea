\name{genperms}
\alias{genperms}
\title{Generate permutations of a phenotype vector}
\usage{
  genperms(pheno, n.perm = 0)
}
\arguments{
  \item{pheno}{a vector to be permuted}

  \item{n.perm}{the number of times to permute}
}
\value{
  a matrix with dimensions length(pheno) x n.perm.
}
\description{
  Given a vector, generate n.perm samples and return a
  matrix with each permutation in each column.
}
\examples{
y <- rbinom(50,2,0.3)
genperms(y,4)
}
\author{
  Chris Wallace <chris.wallace at cimr.cam.ac.uk>
}
\keyword{manip}

