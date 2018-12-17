\name{setParams}
\alias{setParams}
\title{Parameter setting for lasso}
\description{
This functions initializes all required parameters for lasso and tree lasso to defaults.
}
\usage{
setParams(opts)
}
\arguments{
\item{opts}{Either an empty list or list of options with values.}
}
\value{
	\item{opts}{Default parameter settings.}
}
\examples{
opts <- list()
opts <- setParams(opts)
}
\keyword{file}
