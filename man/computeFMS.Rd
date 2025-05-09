% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeFMS.R
\name{computeFMS}
\alias{computeFMS}
\title{Compute Factor Match Score for two models.}
\usage{
computeFMS(Fac1, Fac2, modes)
}
\arguments{
\item{Fac1}{A list of matrices corresponding to found components per mode in model 1.}

\item{Fac2}{A list of matrices corresponding to found components per mode in model 2.}

\item{modes}{List of modes per dataset.}
}
\value{
Vector of FMS scores, one per dataset.
}
\description{
Compute Factor Match Score for two models.
}
\examples{
A = array(rnorm(108*2), c(108, 2))
B = array(rnorm(100*2), c(100, 2))
C = array(rnorm(10*2), c(10, 2))
D = array(rnorm(100*2), c(100, 2))
E = array(rnorm(10*2), c(10, 2))

Fac1 = list(A,B,C,D,E)
Fac2 = Fac1 # identical models for the purposes of demonstration
modes = list(c(1,2,3), c(1,4,5))
FMS_result = computeFMS(Fac1, Fac2, modes) # FMS_result = c(1,1)
}
