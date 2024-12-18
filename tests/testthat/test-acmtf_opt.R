test_that("a solution is found in the two-tensor case", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)

  expect_no_error(acmtf_opt(Z, 1, max_iter=2))
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

  result = acmtf_opt(Z, 2, initialization="nvec")
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

  result = acmtf_opt(Z, 2, initialization="random", max_iter = 2)
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

  results = acmtf_opt(Z, 2, initialization="random", nstart=10, max_iter=2, allOutput=TRUE)
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

  expect_no_error(acmtf_opt(Z,2,initialization="random", nstart=2, max_iter=2, numCores=2))
})

test_that("acmtf_opt works when including a vector Y (1 comp)", {
  numComponents = 1
  I = 108
  J = 100
  K = 10
  L = 1
  A = array(rnorm(I*numComponents), c(I, numComponents))  # shared subject mode
  B = array(rnorm(J*numComponents), c(J, numComponents))  # distinct feature mode of X1
  C = array(rnorm(K*numComponents), c(K, numComponents))  # distinct condition mode of X1
  D = as.matrix(1)
  lambdas = as.matrix(c(1, 1))

  df1 = array(0L, c(I, J, K))
  df2 = array(0L, c(I, L))
  for(i in 1:numComponents){
    df1 = df1 + lambdas[1,i] * reinflateTensor(A[,i], B[,i], C[,i])
    df2 = df2 + lambdas[2,i] * reinflateMatrix(A[,i], D[,i])
  }
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4))
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)

  expect_no_error(acmtf_opt(Z,1,initialization="random", max_iter=2))
})

test_that("acmtf_opt works when including a vector Y (3 comp)", {
  numComponents = 3
  I = 108
  J = 100
  K = 10
  L = 1
  A = array(rnorm(I*numComponents), c(I, numComponents))  # shared subject mode
  B = array(rnorm(J*numComponents), c(J, numComponents))  # distinct feature mode of X1
  C = array(rnorm(K*numComponents), c(K, numComponents))  # distinct condition mode of X1
  D = t(as.matrix(c(1,1,1)))
  lambdas = matrix(c(1,1,1,1,1,1),nrow=2,ncol=numComponents)

  df1 = array(0L, c(I, J, K))
  df2 = array(0L, c(I, L))
  for(i in 1:numComponents){
    df1 = df1 + lambdas[1,i] * reinflateTensor(A[,i], B[,i], C[,i])
    df2 = df2 + lambdas[2,i] * reinflateMatrix(A[,i], D[,i])
  }
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4))
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)

  expect_no_error(acmtf_opt(Z,numComponents,initialization="random", max_iter=2))
})
