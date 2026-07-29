# Set up datasets for (A)CMTF input

Set up datasets for (A)CMTF input

## Usage

``` r
setupCMTFdata(datasets, modes, normalize = TRUE)
```

## Arguments

- datasets:

  List of arrays of datasets. Multi-way and two-way may be combined.

- modes:

  Numbered modes per dataset in a list. Example element 1: 1 2 3 and
  element 2: 1 4 for the X tensor and Y matrix case with a shared
  subject mode.

- normalize:

  Boolean specifying if the datasets should be normalized to Frobenium
  norm 1.

  Note: this function puts zeroes in positions with missing values. The
  indices of missing data are conserved in the output.

## Value

Z, a list with "object" listing the datasets, "sizes" with their size,
"norms" with their norms and "missing" stating the missing data.

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
```
