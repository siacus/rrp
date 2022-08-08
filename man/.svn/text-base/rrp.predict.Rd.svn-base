\name{rrp.predict}
\alias{rrp.predict}
\title{Nonparametric nearest neighbor predictor using RRP dissimilarity}
\description{
This method is a simple nearest neighbor based on the RRP dissimilarity matrix.
It can be used also as a supervised method, i.e. using the outcome variable
in the construction of the RRP dissimilarity matrix and by imputing missing 
values to the test set.
}
\usage{
rrp.predict(x, y, train, test, k = 1)
}
\arguments{
  \item{x}{a \code{XPtr} object}
  \item{y}{numeric vector of train outcomes}
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
a vector of type \code{numeric} with predicted outcomes.
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
\seealso{\code{\link{rrp.dist}}, \code{\link{rrp.class}}}
\examples{
require(MASS)
attach(birthwt)

race <- factor(race, labels= c("white", "black", "other"))
ptd <- factor(ptl>0)
ftv <- factor(ftv)
levels(ftv)[-(1:2)] <- "2+"
table(ftv)

bwt <- data.frame(bwt, age, lwt, race, smoke=(smoke>0), ptd, 
ht=(ht>0), ui = (ui>0), ftv)

detach()
rm(race, ptd, ftv)

set.seed(123)
n <- dim(bwt)[1]
test <- sample(1:n, 15)
train <- (1:n)[-test]

D <- rrp.dist(bwt[,-1])
true.wht <- bwt$bwt[test]
pred.wht <- rrp.predict(D, bwt$bwt[train], train, test)
mean(pred.wht-true.wht)
sd(pred.wht-true.wht)

mod <- lm(bwt ~ ., data=bwt[train,])
pred.wht <- predict(mod, newdata = bwt[test,])
mean(pred.wht-true.wht)
sd(pred.wht-true.wht)
}
\keyword{tree}
