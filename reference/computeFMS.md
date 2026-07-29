# Compute Factor Match Score for two models.

Compute Factor Match Score for two models.

## Usage

``` r
computeFMS(Fac1, Fac2, modes)
```

## Arguments

- Fac1:

  A list of matrices corresponding to found components per mode in model
  1.

- Fac2:

  A list of matrices corresponding to found components per mode in model
  2.

- modes:

  List of modes per dataset.

## Value

Vector of FMS scores, one per dataset.

## Examples

``` r
A = array(rnorm(108*2), c(108, 2))
B = array(rnorm(100*2), c(100, 2))
C = array(rnorm(10*2), c(10, 2))
D = array(rnorm(100*2), c(100, 2))
E = array(rnorm(10*2), c(10, 2))

Fac1 = list(A,B,C,D,E)
Fac2 = Fac1 # identical models for the purposes of demonstration
modes = list(c(1,2,3), c(1,4,5))
FMS_result = computeFMS(Fac1, Fac2, modes) # FMS_result = c(1,1)
```
