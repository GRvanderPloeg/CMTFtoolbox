# Normalize all vectors in model output Fac object to norm 1.

Normalize all vectors in model output Fac object to norm 1.

## Usage

``` r
normalizeFac(Fac, modes)
```

## Arguments

- Fac:

  List object with all components per mode per item.

- modes:

  List object with modes per dataset (see also
  [`setupCMTFdata()`](https://grvanderploeg.com/CMTFtoolbox/reference/setupCMTFdata.md))

## Value

List object of normalized Fac object, the extracted norms per loading
vector per component, and the norms per dataset per component.

## Examples

``` r
set.seed(123)
A = array(rnorm(108*2), c(108, 2))
B = array(rnorm(100*2), c(100, 2))
C = array(rnorm(10*2), c(10, 2))
D = array(rnorm(100*2), c(100,2))
E = array(rnorm(10*2), c(10,2))
modes = list(c(1,2,3), c(1,4,5))

Fac = list(A, B, C, D, E)
output = normalizeFac(Fac, modes)
```
