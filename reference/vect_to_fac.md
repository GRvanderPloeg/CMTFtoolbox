# Convert vectorized output of (a)cmtf to a Fac list object with all loadings per mode.

Convert vectorized output of (a)cmtf to a Fac list object with all
loadings per mode.

## Usage

``` r
vect_to_fac(vect, Z, sortComponents = FALSE)
```

## Arguments

- vect:

  Vectorized output of (a)cmtf

- Z:

  Original Z input object (see
  [setupCMTFdata](https://grvanderploeg.com/CMTFtoolbox/reference/setupCMTFdata.md)).

- sortComponents:

  Sort the order of the components by variation explained (default
  FALSE).

## Value

Fac: list object with all loadings in all components per mode, ordered
the same way as Z\$modes.

## Examples

``` r
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

result = cmtf_opt(Z, 2, initialization="random", max_iter = 2)
Fac = vect_to_fac(result$par, Z)
```
