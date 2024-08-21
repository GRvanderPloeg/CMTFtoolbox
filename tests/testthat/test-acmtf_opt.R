test_that("a solution is found in the two-tensor case", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)

  expect_no_error(acmtf_opt(Z, 1))
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

  result = acmtf_opt(Z, 2, initialization="nvec")
  expect_equal(result$f, 0, tolerance = 0.1)
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

  result = acmtf_opt(Z, 2, initialization="random", max_iter = 2)
  expect_gt(result$f, 1)
})

# test_that("the lambdas cannot be negative", {
#   set.seed(123)
#   numComponents = 3
#   I = 108
#   J = 100
#   K = 10
#   L = 100
#   A = array(rnorm(I*numComponents), c(I, numComponents))  # shared subject mode
#   B = array(rnorm(J*numComponents), c(J, numComponents))  # distinct feature mode of X1
#   C = array(rnorm(K*numComponents), c(K, numComponents))  # distinct condition mode of X1
#   D = array(rnorm(L*numComponents), c(L, numComponents))  # distinct feature mode of X2
#   lambdas = array(c(1, 1, 1, 0, 0, 1), c(2,3))
#
#   df1 = array(0L, c(I, J, K))
#   df2 = array(0L, c(I, L))
#   for(i in 1:numComponents){
#     df1 = df1 + lambdas[1,i] * reinflateTensor(A[,i], B[,i], C[,i])
#     df2 = df2 + lambdas[2,i] * reinflateMatrix(A[,i], D[,i])
#   }
#   datasets = list(df1, df2)
#   modes = list(c(1,2,3), c(1,4))
#   Z = setupCMTFdata(datasets, modes, normalize=TRUE)
#
#   result = acmtf_opt(Z, 3, initialization="nvec")
#
#   expect_true(all(result$Fac[[5]] >= 0))
# })
