# Create a matrix from a matrix of scores and loadings similar to a component model.

Create a matrix from a matrix of scores and loadings similar to a
component model.

## Usage

``` r
reinflateMatrix(A, B)
```

## Arguments

- A:

  I x N matrix corresponding to scores for N components.

- B:

  J x N matrix corresponding to loadings for N components.

## Value

M, an I x J matrix.

## Examples

``` r
A = rnorm(108)
B = rnorm(100)
M = reinflateMatrix(A,B)
```
