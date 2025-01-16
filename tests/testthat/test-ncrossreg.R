test_that("the function throws no errors with normal use", {
  set.seed(123)
  A = array(rnorm(25*2), c(25, 2))
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

  expect_no_error(ncrossreg(Z, Y, maxNumComponents=2, max_iter=2, nstart=2))
})

test_that("running in parallel works", {
  skip_on_cran()

  set.seed(123)
  A = array(rnorm(25*2), c(25, 2))
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

  expect_no_error(ncrossreg(Z, Y, maxNumComponents=2, max_iter=2, nstart=2, numCores=2))
})
