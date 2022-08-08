\name{newXPtr}
\alias{newXPtr}
\title{Creates a new XPtr object}
\description{
Creates a new \code{XPtr} of size \code{n} and initialize it with the constant
\code{k}. Objects of class \code{XPtr} are just \code{externalptr} with
attributes. They are used instead of classical \code{dist} to avoid copies. 
This makes RRP to work on bigger data sets.
}
\usage{
newXPtr(n, k = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{the size of the object}
  \item{k}{the initializing constant}
}
\details{
You should think about \code{XPtr} as \code{dist} objects, i.e. a representation
of a symmetric matrix with 0's on the diagonals. \code{XPtr} are stored as
lower diagonal matrices. The size \code{n} is the dimension of the corresponding
matrix, say \code{M}, hence \code{M} is intended as a \code{n} * \code{n} 
matrix. The real length of the \code{XPtr} is \code{n*(n-1)/2}.
}
\value{
returns a new \code{XPtr} if there is enough memory.
}
\references{Iacus, S.M., Porro, G. (2006) Random Recursive Partitioning and 
its applications to missing data imputation, classification and 
average treatment effect estimation, \emph{submitted}.}
\author{S.M. Iacus}
\seealso{\code{\link{rrp.dist}}, \code{\link{XPtrToDist}}}
\examples{
a <- newXPtr(10, 1)
(XPtrToDist(a))
as.dist(matrix(1,10,10))
}
\keyword{tree}
