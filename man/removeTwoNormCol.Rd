% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{removeTwoNormCol}
\alias{removeTwoNormCol}
\title{Remove two-norms column-wise from a matrix}
\usage{
removeTwoNormCol(df)
}
\arguments{
\item{df}{Matrix of loadings}
}
\value{
Matrix of loadings where the column-wise 2-norm is 1.
}
\description{
Remove two-norms column-wise from a matrix
}
\examples{
A = array(rnorm(108*4), c(108,4))
Anorm = removeTwoNormCol(A)
}
