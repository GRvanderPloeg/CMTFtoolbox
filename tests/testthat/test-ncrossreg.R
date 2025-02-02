test_that("the function throws no errors with normal use", {
  set.seed(123)
  I = 10
  J = 5
  K = 3
  L = 8
  M = 3
  A = array(rnorm(I*2), c(I, 2))
  B = array(rnorm(J*2), c(J, 2))
  C = array(rnorm(K*2), c(K, 2))
  D = array(rnorm(L*2), c(L, 2))
  E = array(rnorm(M*2), c(M, 2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateTensor(A, D, E)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  Y = matrix(A[,1])

  expect_no_error(ncrossreg(Z, Y, maxNumComponents=2, max_iter=2, nstart=2))
})

test_that("the function throws no errors with K folds", {
  set.seed(123)
  I = 10
  J = 5
  K = 3
  L = 8
  M = 3
  A = array(rnorm(I*2), c(I, 2))
  B = array(rnorm(J*2), c(J, 2))
  C = array(rnorm(K*2), c(K, 2))
  D = array(rnorm(L*2), c(L, 2))
  E = array(rnorm(M*2), c(M, 2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateTensor(A, D, E)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  Y = matrix(A[,1])

  expect_no_error(ncrossreg(Z, Y, maxNumComponents=2, max_iter=2, nstart=2, cvFolds=2))
})

test_that("the function throws an errors if Z is wrong", {
  set.seed(123)
  I = 10
  J = 5
  K = 3
  L = 8
  M = 3
  A = array(rnorm(I*2), c(I, 2))
  B = array(rnorm(J*2), c(J, 2))
  C = array(rnorm(K*2), c(K, 2))
  D = array(rnorm(L*2), c(L, 2))
  E = array(rnorm(M*2), c(M, 2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateTensor(A, D, E)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  Y = matrix(A[,1])

  expect_error(ncrossreg(NULL, Y, maxNumComponents=2, max_iter=2, nstart=2, cvFolds=2))
})

test_that("the function throws an errors if Y is wrong", {
  set.seed(123)
  I = 10
  J = 5
  K = 3
  L = 8
  M = 3
  A = array(rnorm(I*2), c(I, 2))
  B = array(rnorm(J*2), c(J, 2))
  C = array(rnorm(K*2), c(K, 2))
  D = array(rnorm(L*2), c(L, 2))
  E = array(rnorm(M*2), c(M, 2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateTensor(A, D, E)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  Y = matrix(A[1:10,1])

  expect_error(ncrossreg(NULL, Y, maxNumComponents=2, max_iter=2, nstart=2, cvFolds=2))
})

test_that("running in parallel works", {
  skip_on_cran()

  set.seed(123)
  I = 10
  J = 5
  K = 3
  L = 8
  M = 3
  A = array(rnorm(I*2), c(I, 2))
  B = array(rnorm(J*2), c(J, 2))
  C = array(rnorm(K*2), c(K, 2))
  D = array(rnorm(L*2), c(L, 2))
  E = array(rnorm(M*2), c(M, 2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateTensor(A, D, E)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  Y = matrix(A[,1])

  expect_no_error(ncrossreg(Z, Y, maxNumComponents=2, max_iter=2, nstart=2, numCores=2, cvFolds=2))
})
