% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncrossreg.R
\name{ncrossreg}
\alias{ncrossreg}
\title{Cross-validation of ACMTF-R by classical K-fold CV with best-model selection per fold}
\usage{
ncrossreg(
  Z,
  Y,
  maxNumComponents = 5,
  alpha = 1,
  beta = rep(0.001, length(Z$object)),
  epsilon = 1e-08,
  pi = 0.5,
  normalize = TRUE,
  normY = 1,
  method = "CG",
  cg_update = "HS",
  line_search = "MT",
  max_iter = 10000,
  max_fn = 10000,
  abs_tol = 1e-10,
  rel_tol = 1e-10,
  grad_tol = 1e-10,
  nstart = 5,
  numCores = 1,
  cvFolds = 2
)
}
\arguments{
\item{Z}{Combined dataset and mode object as produced by \code{\link[=setupCMTFdata]{setupCMTFdata()}}.}

\item{Y}{Dependent variable (regression part).}

\item{maxNumComponents}{Maximum number of components to investigate (default 5).}

\item{alpha}{Scalar penalizing the components to be norm 1 (default 1).}

\item{beta}{Vector of penalty values for each dataset, penalizing the lambda terms (default 1e-3).}

\item{epsilon}{Scalar value to make it possible to compute the partial derivatives of lambda (default 1e-8).}

\item{pi}{Pi value of the loss function as specified by Van der Ploeg et al., 2025.}

\item{normalize}{Normalize the X blocks to frobenius norm 1 (default TRUE).}

\item{normY}{Normalize Y to a specific value, (default: 1).}

\item{method}{Optimization method to use (default = "CG", the conjugate gradient). See \code{\link[mize:mize]{mize::mize()}} for other options.}

\item{cg_update}{Update method for the conjugate gradient algorithm, see \code{\link[mize:mize]{mize::mize()}} for the options (default="HS", Hestenes-Steifel).}

\item{line_search}{Line search algorithm to use, see \code{\link[mize:mize]{mize::mize()}} for the options (default="MT", More-Thuente).}

\item{max_iter}{Maximum number of iterations.}

\item{max_fn}{Maximum number of function evaluations.}

\item{abs_tol}{Function tolerance criterion for convergence.}

\item{rel_tol}{Relative function tolerance criterion for convergence.}

\item{grad_tol}{Absolute tolerence for the l2-norm of the gradient vector.}

\item{nstart}{Number of replicate initializations per CV fold (default 5).}

\item{numCores}{Number of cores to use for the replicate fits (default 1). Each replicate model is
fit with \code{numCores = 1} (to avoid nested parallelism) but the replicates are run in parallel.}

\item{cvFolds}{Number of folds to use in the cross-validation. For example, if \code{cvFolds}
is 5, then the subjects are deterministically partitioned into 5 groups
(each CV iteration uses 4/5 for training and 1/5 for testing). Default: 2.}
}
\value{
A list with two elements:
- \strong{varExp}: a tibble with the variance–explained (for X and Y) per number of components.
- \strong{RMSE}: a tibble with the RMSE (computed over the unified CV prediction vector) per number of components.
- \strong{FMS}: a tibble with Factor Match Score per data block per number of components.
}
\description{
This function runs ACMTF-R with cross-validation. A deterministic K–fold partition
is used: the subjects are split in order into \code{cvFolds} groups. For each fold the
training set consists of the other folds and the test set is the current fold.
}
\details{
For each fold and for each number of components, \emph{nstart} models are fitted—each
with a single initialization (i.e. using \code{nstart = 1} and \code{numCores = 1} in \code{acmtfr_opt}).
The large parallel loop iterates over all combinations of (numComponents, fold, replicate).
After the loop, for each (numComponents, fold) the best replicate is selected (using \code{model$f})
and used to predict on the test set. The predictions for all folds are then combined into
a unified prediction vector for RMSE calculation.
}
\examples{
set.seed(123)
A = array(rnorm(25 * 2), c(25, 2))
B = array(rnorm(100 * 2), c(100, 2))
C = array(rnorm(10 * 2), c(10, 2))
D = array(rnorm(100 * 2), c(100, 2))
E = array(rnorm(10 * 2), c(10, 2))

df1 = reinflateTensor(A, B, C)
df2 = reinflateTensor(A, D, E)
datasets = list(df1, df2)
modes = list(c(1, 2, 3), c(1, 4, 5))
Z = setupCMTFdata(datasets, modes)
Y = matrix(A[, 1])

# For classical 5-fold CV (deterministic splits) with best-model selection:
result = ncrossreg(Z, Y, maxNumComponents = 2, max_iter = 2, nstart = 2, cvFolds = 5)
}
