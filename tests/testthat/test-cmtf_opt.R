test_that("a solution is found in the two-tensor case", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)

  expect_no_error(cmtf_opt(Z, 1))
})

test_that("the objective is close to zero if the correct solution is found", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)

  result = cmtf_opt(Z, 2, initialization="nvec")
  expect_equal(result$f, 0, tolerance = 1e-4)
})

test_that("the objective is very high if an incorrect solution is found", {
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)

  result = cmtf_opt(Z, 2, initialization="random", max_iter = 2)
  expect_gt(result$f, 1)
})
