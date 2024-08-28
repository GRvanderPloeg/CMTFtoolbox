test_that("investigateFMS runs without error in acmtf mode", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))

  expect_no_error(investigateFMS(datasets, modes, 1, model="acmtf", numFolds=2, rel_tol=1e-5, abs_tol=1e-5))
})

test_that("investigateFMS runs without error when jack-knifing", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))

  expect_no_error(investigateFMS(datasets, modes, 1, model="acmtf", numFolds=2, jackKnife=TRUE, rel_tol=1e-5, abs_tol=1e-5))
})

test_that("investigateFMS runs without error in cmtf mode", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))

  expect_no_error(investigateFMS(datasets, modes, 1, model="cmtf", numFolds=2, rel_tol=1e-5, abs_tol=1e-5))
})

test_that("specifying the wrong type of model results in an error", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))

  expect_error(investigateFMS(datasets, modes, 1, model="test", numFolds=2, rel_tol=1e-5, abs_tol=1e-5))
})

test_that("running in parallel works", {
  skip_on_cran()

  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))

  expect_no_error(investigateFMS(datasets, modes, 1, model="acmtf", numFolds=2, numCores=2, rel_tol=1e-5, abs_tol=1e-5))
})
