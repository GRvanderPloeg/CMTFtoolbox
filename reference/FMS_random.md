# Compute Factor Match Score for two models.

Compute Factor Match Score for two models.

## Usage

``` r
FMS_random(Fac1, Fac2)
```

## Arguments

- Fac1:

  A list of matrices corresponding to found components per mode in model
  1.

- Fac2:

  A list of matrices corresponding to found components per mode in model
  2.

## Value

Scalar of FMS value

## Examples

``` r
set.seed(123)

I = 10
J = 5
K = 3
df = array(rnorm(I*J*K), c(I,J,K))
datasets = list(df, df)
modes = list(c(1,2,3), c(1,4,5))
Z = setupCMTFdata(datasets, modes)

model1 = acmtf_opt(Z, 1)

Fac1 = model1$Fac[1:3]
Fac2 = Fac1 # identical models for the purposes of demonstration
result = FMS_random(Fac1, Fac2) # [1] 1
```
