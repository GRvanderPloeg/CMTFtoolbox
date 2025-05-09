% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initializeCMTF.R
\name{initializeCMTF}
\alias{initializeCMTF}
\title{Initialize input vectors for the CMTF algorithm}
\usage{
initializeCMTF(
  Z,
  numComponents,
  initialization = "random",
  output = "Fac",
  Y = NULL
)
}
\arguments{
\item{Z}{List object as generated by \code{\link[=setupCMTFdata]{setupCMTFdata()}}.}

\item{numComponents}{Integer stating the number of desired components for the CMTF model.}

\item{initialization}{Initialization method, either "random" or "nvec" (default "random"). Random will initialize random input vectors. Nvec will initialize vectors according to an singular value decomposition of the (matricized, if needed) concatenated datasets per mode.}

\item{output}{How to return output: as a "Fac" object (default) or vectorized ("vect").}

\item{Y}{Used as dependent variable when initialization is set to "npls". Not used by default.}
}
\value{
List or vector of initialized input vectors per mode.
}
\description{
Initialize input vectors for the CMTF algorithm
}
\examples{
set.seed(123)
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

init = initializeCMTF(Z, 1)
}
