test_that("the function works normally", {
  set.seed(123)

  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*2), c(100, 2))
  C = array(rnorm(10*2), c(10, 2))
  D = array(rnorm(100*2), c(100,2))
  E = array(rnorm(10*2), c(10,2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateTensor(A, D, E)
  Y = A[,1]
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)

  expect_no_error(cv_degeneracy(Z, Y, numComponents=1:2, pis=c(0.50, 0.75), rel_tol=1e-4, abs_tol=1e-4, nstart=1, numCores=1))
})
