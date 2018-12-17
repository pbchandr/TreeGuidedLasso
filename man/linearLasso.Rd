\name{linearLasso}
\alias{linearLasso}
\title{Linear Regerssion with lasso}
\description{
This functions perfroms regression with lasso constarints.
}
\usage{
linearLasso(A, y, z, opts)
}
\arguments{
\item{A}{data matrix with samples as rows and features as columns.}
\item{y}{Continuous response values for each sample.}
\item{z}{Regularization parameter. This defines the number of features to be selected.}
\item{opts}{Setting the optional paramters. All parameters in opts need not be defined as it will be set to defaults.}
}
\value{
	\item{ll_out}{Output list containing three elements: beta-weights of each feature, valueL-eigen value, funVal-function value.}
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
z <- 0.01
ll_out <- linearLasso(A, y, z, opts)
}
\keyword{file}
