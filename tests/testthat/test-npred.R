test_that("the function throws no errors with one new sample", {
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

  model = acmtfr_opt(Ztrain,Ytrain,2,initialization="random",pi=0, nstart=1, max_iter=10)
  expect_no_error(npred(model, Xtest, Ztrain))
})

test_that("the function throws no errors with several new samples", {
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
  i = c(1,2)
  Xtest = lapply(Z$object, function(x){x@data[i,,]})
  Ytest = Y[i]

  Xtrain = lapply(Z$object, function(x){x@data[-i,,]})
  Ytrain = Y[-i]
  Ztrain = setupCMTFdata(Xtrain, Z$modes)

  model = acmtfr_opt(Ztrain,Ytrain,2,initialization="random",pi=0, nstart=1, max_iter=10)
  expect_no_error(npred(model, Xtest, Ztrain))
})

test_that("vectX must be the same size as vectZ", {
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
  Xtest = lapply(Z$object, function(x){x@data[i,1:10,]})
  Ytest = Y[i]

  Xtrain = lapply(Z$object, function(x){x@data[-i,,]})
  Ytrain = Y[-i]
  Ztrain = setupCMTFdata(Xtrain, Z$modes)

  model = acmtfr_opt(Ztrain,Ytrain,2,initialization="random",pi=0, nstart=1, max_iter=10)
  expect_error(npred(model, Xtest, Ztrain))
})

test_that("missing values in Xnew are ignored for the prediction", {
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
  Xtest[[1]][1:10,] = NA
  Ytest = Y[i]

  Xtrain = lapply(Z$object, function(x){x@data[-i,,]})
  Ytrain = Y[-i]
  Ztrain = setupCMTFdata(Xtrain, Z$modes)

  model = acmtfr_opt(Ztrain,Ytrain,2,initialization="random",pi=0, nstart=1, max_iter=10)
  expect_error(npred(model, Xtest, Ztrain))
})
