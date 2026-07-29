# Advanced coupled matrix and tensor factorizations

Advanced coupled matrix and tensor factorizations

## Usage

``` r
acmtf_opt(
  Z,
  numComponents,
  initialization = "random",
  alpha = 1,
  beta = rep(0.001, length(Z$object)),
  epsilon = 1e-08,
  method = "CG",
  cg_update = "HS",
  line_search = "MT",
  max_iter = 10000,
  max_fn = 10000,
  abs_tol = 1e-10,
  rel_tol = 1e-10,
  grad_tol = 1e-10,
  nstart = 1,
  numCores = 1,
  sortComponents = TRUE,
  allOutput = FALSE
)
```

## Arguments

- Z:

  Combined dataset and mode object as produced by
  [`setupCMTFdata()`](https://grvanderploeg.com/CMTFtoolbox/reference/setupCMTFdata.md).

- numComponents:

  Number of components

- initialization:

  Initialization, either "random" (default) or "nvec" for numComponents
  components of the concatenated data using svd.

- alpha:

  Scalar penalizing the components to be norm 1 (default 1).

- beta:

  Vector of penalty values for each dataset, penalizing the lambda terms
  (default 1e-3).

- epsilon:

  Scalar value to make it possible to compute the partial derivatives of
  lambda (default 1e-8).

- method:

  Optimization method to use (default = "CG", the conjugate gradient).
  See [`mize::mize()`](https://rdrr.io/pkg/mize/man/mize.html) for other
  options.

- cg_update:

  Update method for the conjugate gradient algorithm, see
  [`mize::mize()`](https://rdrr.io/pkg/mize/man/mize.html) for the
  options (default="HS", Hestenes-Steifel).

- line_search:

  Line search algorithm to use, see
  [`mize::mize()`](https://rdrr.io/pkg/mize/man/mize.html) for the
  options (default="MT", More-Thuente).

- max_iter:

  Maximum number of iterations.

- max_fn:

  Maximum number of function evaluations.

- abs_tol:

  Function tolerance criterion for convergence.

- rel_tol:

  Relative function tolerance criterion for convergence.

- grad_tol:

  Absolute tolerence for the l2-norm of the gradient vector.

- nstart:

  Number of models to produce (default 1). If set higher than one, the
  package will return the best fitted model.

- numCores:

  Number of cores to use (default 1). If set higher than one, the
  package will attempt to run in parallel.

- sortComponents:

  Sort the components in the output by descending order of variation
  explained.

- allOutput:

  Return all created models. Ignored if nstart=1.

## Value

List object, similar to
[`mize::mize()`](https://rdrr.io/pkg/mize/man/mize.html) output.
Includes a Fac object of the model, which is a list of components per
mode. Also includes an init object giving the initialized input vectors.

## Examples

``` r
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

# specific setting to reduce runtime for CRAN
model = acmtf_opt(Z, 1, rel_tol=1e-5, abs_tol=1e-5)
```
