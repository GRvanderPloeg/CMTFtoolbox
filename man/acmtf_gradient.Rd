% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acmtf_gradient.R
\name{acmtf_gradient}
\alias{acmtf_gradient}
\title{Calculate gradient of ACMTF model.}
\usage{
acmtf_gradient(
  x,
  Z,
  alpha = 1,
  beta = rep(0.001, length(Z$object)),
  epsilon = 1e-08
)
}
\arguments{
\item{x}{Vectorized parameters of the CMTF model.}

\item{Z}{Z object as generated by \code{\link[=setupCMTFdata]{setupCMTFdata()}}.}

\item{alpha}{Alpha value of the loss function as specified by Acar et al., 2014}

\item{beta}{Beta value of the loss function as specified by Acar et al., 2014}

\item{epsilon}{Epsilon value of the loss function as specified by Acar et al., 2014}
}
\value{
Vectorized gradient of the ACMTF model.
}
\description{
Calculate gradient of ACMTF model.
}
\examples{
A = array(rnorm(108*2), c(108, 2))
B = array(rnorm(100*2), c(100, 2))
C = array(rnorm(10*2), c(10, 2))
D = array(rnorm(100*2), c(100,2))
E = array(rnorm(10*2), c(10,2))

df1 = reinflateTensor(A, B, C)
df2 = reinflateTensor(A, D, E)
datasets = list(df1, df2)
modes = list(c(1,2,3), c(1,4,5))
Z = setupCMTFdata(datasets, modes, normalize=FALSE)

init = initializeACMTF(Z, 2, output="vect")
g = acmtf_gradient(init, Z)
}
