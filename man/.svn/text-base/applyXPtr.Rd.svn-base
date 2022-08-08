\name{applyXPtr}
\alias{applyXPtr}
\title{Apply facility for objects of class XPtr}
\description{
This tool function allow to apply a function along some index and on a subset
of elements of the \code{XPtr} object. Is the counterpart of \code{apply}
for objects fo class \code{matrix}.
}
\usage{
applyXPtr(d, idx, sub, f)
}
\arguments{
  \item{d}{an object of class \code{XPtr}}
  \item{idx}{a vector of integers, the indexes of the apply}
  \item{sub}{a vector of integers, the indexes of the subest}
  \item{f}{a function to be evaluated in the \code{.GlobalEnv}}
}
\details{
The object \code{d} is like a \code{dist} object, which is one-dimensional but
usually associated to a symmetrix matrix \code{M} with 0's on the diagonal.
This function allows to apply a function \code{f} on \code{d} as if it was
applied this way on the matrix \code{M}: \code{apply(M[idx,sub],1,f)}.
}
\value{
  \item{val}{an invisible object of type \code{list}}
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
a <- newXPtr(10,1)
addXPtr(a, list(c(1,3,4,9), c(2,5), c(1,7,9)), c(-1, 5, 10))
f <- function(x) sum(x)
idx <- 1:5
sub <- 7:10
applyXPtr(a, idx, sub, f) -> l
(XPtrToDist(a))
l
}
\keyword{tree}
