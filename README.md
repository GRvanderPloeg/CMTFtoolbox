
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CMTFtoolbox

<!-- badges: start -->

[![codecov](https://codecov.io/gh/GRvanderPloeg/CMTFtoolbox/graph/badge.svg?token=Y8XWFAV0IC)](https://codecov.io/gh/GRvanderPloeg/CMTFtoolbox)
[![R-CMD-check](https://github.com/GRvanderPloeg/CMTFtoolbox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GRvanderPloeg/CMTFtoolbox/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The `CMTFtoolbox` package provides R users with two data fusion methods
that have previously been presented in the MATLAB sphere.

- `cmtf_opt`: Coupled Matrix and Tensor Factorization (CMTF) as
  described in [Acar et al., 2011](https://arxiv.org/abs/1105.3422).
- `acmtf_opt`: Advanced Coupled Matrix and Tensor Factorization (ACMTF)
  as described in [Acar et al.,
  2014](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-239).

Both of these methods were implemented using the all-at-once
optimization approaches as described in the papers above. This
implementation was achieved using the S4 Tensor object from `rTensor`
and the various conjugate gradient approaches from `mize`.

## Installation

You can install the development version of `CMTFtoolbox` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GRvanderPloeg/CMTFtoolbox")
```

## Citation

Pending

## Usage

``` r
library(CMTFtoolbox)

set.seed(123)
numComponents = 3
I = 108
J = 100
K = 10
L = 100
A = array(rnorm(I*numComponents), c(I, numComponents))  # shared subject mode
B = array(rnorm(J*numComponents), c(J, numComponents))  # distinct feature mode of X1
C = array(rnorm(K*numComponents), c(K, numComponents))  # distinct condition mode of X1
D = array(rnorm(L*numComponents), c(L, numComponents))  # distinct feature mode of X2
lambdas = array(c(1, 1, 1, 0, 0, 1), c(2,3))

df1 = array(0L, c(I, J, K))
df2 = array(0L, c(I, L))
for(i in 1:numComponents){
  df1 = df1 + lambdas[1,i] * reinflateTensor(A[,i], B[,i], C[,i])
  df2 = df2 + lambdas[2,i] * reinflateMatrix(A[,i], D[,i])
}
datasets = list(df1, df2)
modes = list(c(1,2,3), c(1,4))
Z = setupCMTFdata(datasets, modes, normalize=TRUE)

cmtf_model = cmtf_opt(Z, 3)
acmtf_model = acmtf_opt(Z, 3)
```
