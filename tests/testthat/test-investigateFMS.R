test_that("investigateFMS runs without error in acmtf mode", {
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

  expect_no_error(investigateFMS(datasets, modes, 1, model="acmtf", numFolds=2, max_iter=2))
})

test_that("investigateFMS runs without error when jack-knifing", {
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

  expect_no_error(investigateFMS(datasets, modes, 1, model="acmtf", numFolds=2, jackKnife=TRUE, max_iter=2))
})

test_that("investigateFMS runs without error in cmtf mode", {
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

  expect_no_error(investigateFMS(datasets, modes, 1, model="cmtf", numFolds=2, max_iter=2))
})

test_that("specifying the wrong type of model results in an error", {
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

  expect_error(investigateFMS(datasets, modes, 1, model="test", numFolds=2, max_iter=2))
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
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)

  expect_no_error(investigateFMS(datasets, modes, 1, model="acmtf", numFolds=2, numCores=2, max_iter=2))
})
