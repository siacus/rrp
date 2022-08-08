\name{rrp.dist}
\alias{rrp.dist}
\title{RRP dissimilarity matrix}
\description{
Main piece of the RRP algorithm which generates the RRP dissimilarity
matrix.
}
\usage{
rrp.dist(X, treated = NULL, msplit = 10, Rep = 250, cut.in = 15, 
 check.bal = FALSE, plot = FALSE, asdist = FALSE, verbose = 0)
}
\arguments{
  \item{X}{a \code{data.frame} object}
  \item{treated}{optional class indicator variable}
  \item{msplit}{minimum split parameter in the \code{rpart} algorithm}
  \item{Rep}{number of RRP replications}
  \item{cut.in}{number of breaks in which to cut continuous variables}
  \item{check.bal}{indicator function. If \code{TRUE} balance check using
  hyper-rectangles will be used inside leaves}
  \item{plot}{wheter to plot the porximity matrix as image}
  \item{asdist}{if \code{TRUE} returns an object of class \code{dist}}
  \item{verbose}{if greater than 1 some information is printed}
}
\details{
This algorithm allows for missing data in \code{X}. From version 1.6 of the
package the RRP matrix is stored as an external pointer to avoid duplications.
This allow to work on bigger datasets. 
}
\value{an object of class \code{externalptr, XPtr} or \code{dist}
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
\seealso{\code{\link{rpart}}}
\keyword{tree}
