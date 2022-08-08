\name{setDist}
\alias{setDist}
\title{Sets the elements of a dist object to a constant vector-wise}
\description{
Sets the elements of a dist object to a constant vector-wise. The idea is that 
a symmetrix matrix \code{M} can be represented as a triangular object \code{d}. 
With this function you set \code{M[x,x] = k} but instead of working on 
\code{M} this function works on its \code{dist} representation \code{d}.
}
\usage{
setDist(d, x, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{a dist object, upper diagonal, representing a symmetric matrix}
  \item{x}{a list of vectors of indices}
  \item{k}{a vector of constants to add, arguments are eventually recycled 
  silently}
}
\value{
  returns the modified dist object
}
\author{S. M. Iacus}
\seealso{\code{\link{dist}}}
\examples{
M <- matrix(0,5,5)
d <- as.dist(M)
x <- list(1:3,4:5)

setDist(d,x,c(1,2))
# which is the equivalent of (apart for the diagonal elements)
M[1:3,1:3] <- 1
M[4:5,4:5] <- 2


# Indeed, we have
d
as.dist(M)
}
\keyword{misc}
