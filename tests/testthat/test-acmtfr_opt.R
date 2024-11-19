test_that("a solution is found in the two-tensor case and Y", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  Y = matrix(rnorm(108), nrow=108, ncol=1)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)

  expect_no_error(acmtfr_opt(Z, Y, 1, max_iter=2))
})

test_that("the objective is close to zero if the correct solution is found", {
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
  Y = A[,1]

  result = acmtfr_opt(Z, Y, 2, initialization="nvec")
  expect_equal(result$f, 0, tolerance = 0.1)
})

test_that("the objective is very high if an incorrect solution is found", {
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
  Y = matrix(rnorm(108), nrow=108, ncol=1)

  result = acmtfr_opt(Z, Y, 2, initialization="random", max_iter = 2)
  expect_gt(result$f, 1)
})

test_that("allOutput=TRUE gives a list of expected length", {
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
  Y = matrix(rnorm(108), nrow=108, ncol=1)

  results = acmtfr_opt(Z, Y, 2, initialization="random", nstart=10, max_iter=2, allOutput=TRUE)
  expect_equal(length(results), 10)
})

test_that("running in parallel works", {
  skip_on_cran()

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
  Y = matrix(rnorm(108), nrow=108, ncol=1)

  expect_no_error(acmtfr_opt(Z,Y,2,initialization="random", nstart=2, max_iter=2, numCores=2))
})
