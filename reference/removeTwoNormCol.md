# Remove two-norms column-wise from a matrix

Remove two-norms column-wise from a matrix

## Usage

``` r
removeTwoNormCol(df)
```

## Arguments

- df:

  Matrix of loadings

## Value

Matrix of loadings where the column-wise 2-norm is 1.

## Examples

``` r
A = array(rnorm(108*4), c(108,4))
Anorm = removeTwoNormCol(A)
```
