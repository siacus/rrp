\name{mulXPtr}
\alias{mulXPtr}
\title{Multiplies the elements of a XPtr object by a constant vector-wise}
\description{
Multiplies the elements of a \code{XPtr} object by a constant vector-wise. The 
idea is  that a symmetrix matrix \code{M} can be represented as a triangular 
object \code{d}. With this function you set \code{M[x,x] = M[x,x] * k} but 
instead of working on \code{M} this function works on its \code{XPtr} 
representation \code{d}.
}
\usage{
mulXPtr(d, x, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{a \code{XPtr} object, lower diagonal, representing a symmetric 
  matrix}
  \item{x}{a list of vectors of indices}
  \item{k}{a vector of multiplicative constants, arguments are eventually recycled 
  silently}
}
\value{
  returns the modified \code{XPtr} object
}
\author{Stefano M. Iacus}
\seealso{\code{\link{dist}}}
\examples{
M <- matrix(2,5,5)
d <- newXPtr(5,2)
x <- list(1:3, 4:5)

mulXPtr(d,x,c(-1,2))

# which is the equivalent of (apart for the diagonal elements)
M[1:3,1:3] <- M[1:3,1:3] * (-1)
M[4:5,4:5] <- M[4:5,4:5] * 2

# Indeed, we have
(XPtrToDist(d))
as.dist(M)
}
\keyword{misc}
