\name{rrp.class}
\alias{rrp.class}
\title{Nonparametric nearest neighbor classification using RRP dissimilarity}
\description{
This method is a simple nearest neighbor based on the RRP dissimilarity matrix.
It can be used also as a supervised method, i.e. using the class variable
in the construction of the RRP dissimilarity matrix and by imputing missing 
values to the test set.
}
\usage{
rrp.class(x, cl, train, test, k = 1)
}
\arguments{
  \item{x}{a \code{XPtr} object}
  \item{cl}{class vector, coerced to be of type \code{factor}}
  \item{train}{the vecotr of indexes of the training set}
  \item{test}{the vector of training indexes of test set}
  \item{k}{number of nearest to consider}
}
\details{
From version 1.6 of the package the RRP matrix is stored as an external pointer 
to avoid duplications. This allow to work on bigger datasets.
Hence this function no longer accepts \code{dist} objects.
}
\value{
a vector of type \code{factor} with predicted classes.
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
\seealso{\code{\link{rrp.dist}}}
\examples{
data(iris)
X <- iris[,-5]
n <- dim(X)[1]

set.seed(123)
test <- sample(1:n, 10)
train <- (1:n)[-test]

## unsupervised
D <- rrp.dist(X)
pred <- rrp.class(D, iris[train,5], train, test)
table(pred, iris[test, 5])

# supervised
X <- iris
X[test,5] <- NA
D <- rrp.dist(X)
pred <- rrp.class(D, iris[train,5], train, test)
table(pred, iris[test, 5])
}
\keyword{classif}
