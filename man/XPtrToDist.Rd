\name{XPtrToDist}
\alias{XPtrToDist}
\title{Coercion from an XPtr object into a dist one}
\description{
Coercies  an \code{XPtr} object into a \code{dist} object which can be
used in other algorithms.
}
\usage{
XPtrToDist(d)
}
\arguments{
  \item{d}{an \code{XPtr} object}
}
\details{
This routine does not duplicate but just makes available as a \code{dist}
 the  data contained in \code{d}
}
\value{returns an invisible object of class \code{dist}
}
\references{
Iacus, S.M., Porro, G. (2009) Random Recursive Partitioning: a matching 
method for the estimation of the average treatment effect, 
\emph{Journal of Applied Econometrics}, 24, 163-185. 

Iacus, S.M., Porro, G. (2007) Missing data imputation, matching and other 
applications of random recursive partitioning, 
\emph{Computational Statistics and Data Analysis}, 52, 2, 773-789.
}
\author{S.M. Iacus}
\seealso{\code{\link{newXPtr}}}
\examples{
a <- newXPtr(10,5)
(XPtrToDist(a))
as.dist(matrix(5,10,10))
}
\keyword{tree}
