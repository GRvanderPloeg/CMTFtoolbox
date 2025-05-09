% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmtf_gradient.R
\name{cmtf_gradient}
\alias{cmtf_gradient}
\title{Calculate gradient of CMTF model.}
\usage{
cmtf_gradient(x, Z)
}
\arguments{
\item{x}{Vectorized parameters of the CMTF model.}

\item{Z}{Z object as generated by \code{\link[=setupCMTFdata]{setupCMTFdata()}}.}
}
\value{
Vectorized gradient of the CMTF model.
}
\description{
Calculate gradient of CMTF model.
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

init = initializeCMTF(Z, 2, output="vect")
g = cmtf_gradient(init, Z)
}
