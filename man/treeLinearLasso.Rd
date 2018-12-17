\name{treeLinearLasso}
\alias{treeLinearLasso}
\title{Tree Guided Linear Regerssion}
\description{
This functions perfroms regression with lasso constarints on tree structured data.
}
\usage{
treeLinearLasso(A, y, z, opts)
}
\arguments{
\item{A}{data matrix with samples as rows and features as columns.}
\item{y}{Continuous response values for each sample.}
\item{z}{Regularization parameter. This defines the number of features to be selected.}
\item{opts}{Setting the optional paramters. All parameters in opts need not be defined except 'ind'
parameter. 'ind' defines the tree structure of the data.}
}
\value{
	\item{tll_out}{Output list containing three elements: beta-weights of each feature, valueL-eigen value, funVal-function value.}
}
\examples{
M <- 50
N <- 100

A <- matrix(rnorm(N*M), M, N)
y <- A\%*\%runif(N) + runif(M)*0.01


opts <- list()
# Starting point
opts[["init"]] <- 2 # Starting from a zero point
# Termination Criteria
opts[["tFlag"]] <- 6
opts[["maxIter"]] <- 100
# Regularization
opts[["rFlag"]] <- 1 # the input parameter 'rho' is a ratio in (0, 1)
opts[["ind"]]= t(rbind( t(c(-1, -1, 1)), t(c(1, 20, sqrt(20))), t(c(21, 40, sqrt(20))), t(c(41, 50, sqrt(10))), t(c(51, 70, sqrt(20))), t(c(71,100, sqrt(30))), t(c(1, 50, sqrt(50))), t(c(51, 100, sqrt(50)))))

z <- 0.01
tll_out <- treeLinearLasso(A, y, z, opts)

}
\keyword{file}
