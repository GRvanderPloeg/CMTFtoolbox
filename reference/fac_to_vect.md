# Vectorize Fac object

Vectorize Fac object

## Usage

``` r
fac_to_vect(Fac)
```

## Arguments

- Fac:

  Fac object from CMTF and ACMTF

## Value

Vectorized Fac object

## Examples

``` r
set.seed(123)
A = array(rnorm(108*2), c(108, 2))
B = array(rnorm(100*2), c(100, 2))
C = array(rnorm(10*2), c(10, 2))
D = array(rnorm(100*2), c(100,2))
E = array(rnorm(10*2), c(10,2))
Fac = list(A, B, C, D, E)
v = fac_to_vect(Fac)
```
