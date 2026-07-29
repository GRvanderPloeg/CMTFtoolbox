# Predict Y for new data by projecting the data onto the latent space defined by an ACMTF-R model.

Predict Y for new data by projecting the data onto the latent space
defined by an ACMTF-R model.

## Usage

``` r
npred(model, newX, Z, sharedMode = 1)
```

## Arguments

- model:

  ACMTF-R model

- newX:

  List object of new data, where each element corresponds to a block

- Z:

  Original input data used for the model

- sharedMode:

  Shared mode between the blocks (default 1).

## Value

Ypred: the predicted value of Y for the new data

## Examples

``` r
set.seed(123)
A = array(rnorm(108*2), c(108, 2))
B = array(rnorm(100*2), c(100, 2))
C = array(rnorm(10*2), c(10, 2))
D = array(rnorm(100*2), c(100, 2))
E = array(rnorm(10*2), c(10, 2))

df1 = reinflateTensor(A, B, C)
df2 = reinflateTensor(A, D, E)
datasets = list(df1, df2)
modes = list(c(1,2,3), c(1,4,5))
Z = setupCMTFdata(datasets, modes)
Y = matrix(A[,1])
# Remove a sample and define
i = 1
Xtest = lapply(Z$object, function(x){x@data[i,,]})
Ytest = Y[i]
Xtrain = lapply(Z$object, function(x){x@data[-i,,]})
Ytrain = Y[-i]
Ztrain = setupCMTFdata(Xtrain, Z$modes)
model = acmtfr_opt(Ztrain,Ytrain,1,initialization="random",pi=1, nstart=1, max_iter=10)
Ypred = npred(model, Xtest, Ztrain, sharedMode=1)
```
