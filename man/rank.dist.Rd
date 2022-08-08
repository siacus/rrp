\name{rank.dist}
\alias{rank.dist}
\title{rank-based dissimilarity matrix}
\description{
Algorithm which generates the rank based proximity
matrix.
}
\usage{
rank.dist(X, msplit = 2, cut.in = 0, 
   thr=0.75, weights, asdist=FALSE, verbose = 0)
}
\arguments{
  \item{X}{a \code{data.frame} object}
  \item{msplit}{minimum split parameter}
  \item{cut.in}{number of breaks in which to cut continuous variables}
  \item{thr}{a threshold to keep only the proximities greater than the threshold}
  \item{weights}{a vector of weights}
  \item{asdist}{if \code{TRUE} returns an object of class \code{dist} 
  otherwise a list of the same length of the data frame}
  \item{verbose}{if greater than 1 some information is printed}
}
\details{
This algorithm allows for missing data in \code{X}. 
}
\value{an object of class \code{list} or \code{dist}
}
\author{S.M. Iacus}
\seealso{\code{\link{rrp.dist}}}
\keyword{tree}
