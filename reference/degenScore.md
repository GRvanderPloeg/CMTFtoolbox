# degenScore

Computes the maximum absolute off-diagonal Tucker congruence coefficient
between subject-mode components in an ACMTF model. This metric serves as
a diagnostic tool to detect potential degeneracy in the subject-mode
loadings.

## Usage

``` r
degenScore(A)
```

## Arguments

- A:

  A numeric matrix of subject-mode loadings (dimensions: subjects x
  components).

## Value

A numeric scalar representing the maximum absolute off-diagonal Tucker
congruence between components.

## Details

A high `degenScore` (e.g., \> 0.85) indicates that two or more
components in the subject mode are highly similar, suggesting a possible
degeneracy or lack of uniqueness. A low value (e.g., \< 0.3) indicates
well-separated components.

## Examples

``` r
# Example: Compute degenScore for a random loading matrix
A <- matrix(rnorm(100), nrow = 10, ncol = 10)
degenScore(A)
#> [1] 0.7334973
```
