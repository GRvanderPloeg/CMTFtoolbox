test_that("the correct mode 1 size is returned", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  result = initializeACMTF(Z, 1, initialization="random")
  expect_equal(nrow(result[[1]]), I)
})

test_that("max(modes)+1 number of initialized components are returned", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  result = initializeACMTF(Z, 1, initialization="random")
  expect_equal(length(result), 6)
})

test_that("the correct mode 1 components are found using nvecs", {
  A = array(rnorm(108))
  B = array(rnorm(100*2), c(100, 2))
  C = array(rnorm(10*2), c(10, 2))

  df1 = array(tcrossprod(A, multiway::krprod(as.matrix(C[,1]), as.matrix(B[,1]))), c(108,100,10))
  df2 = array(tcrossprod(A, multiway::krprod(as.matrix(C[,2]), as.matrix(B[,2]))), c(108,100,10))
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  result = initializeACMTF(Z, 1, initialization="nvec")
  expect_equal(abs(cor(result[[1]], A))[1,1], 1, tolerance=0.01)
})

test_that("the correct mode 2 components are found using nvecs", {
  A = array(rnorm(108))
  B = array(rnorm(100))
  C = array(rnorm(10*2), c(10, 2))

  df1 = array(tcrossprod(A, multiway::krprod(as.matrix(C[,1]), as.matrix(B))), c(108,100,10))
  df2 = array(tcrossprod(A, multiway::krprod(as.matrix(C[,2]), as.matrix(B))), c(108,100,10))
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  result = initializeACMTF(Z, 1, initialization="nvec")
  expect_equal(abs(cor(result[[2]], B))[1,1], 1, tolerance=0.01)
})

test_that("the correct mode 3 components are found using nvecs", {
  A = array(rnorm(108))
  B = array(rnorm(100))
  C = array(rnorm(10))

  df1 = array(tcrossprod(A, multiway::krprod(as.matrix(C), as.matrix(B))), c(108,100,10))
  df2 = array(tcrossprod(A, multiway::krprod(as.matrix(C), as.matrix(B))), c(108,100,10))
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  result = initializeACMTF(Z, 1, initialization="nvec")
  expect_equal(abs(cor(result[[3]], C))[1,1], 1, tolerance=0.01)
})

test_that("the correct number of lambda values are produced", {
  A = array(rnorm(108))
  B = array(rnorm(100))
  C = array(rnorm(10))

  df1 = array(tcrossprod(A, multiway::krprod(as.matrix(C), as.matrix(B))), c(108,100,10))
  df2 = array(tcrossprod(A, multiway::krprod(as.matrix(C), as.matrix(B))), c(108,100,10))
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  result = initializeACMTF(Z, 2, initialization="random")
  expect_equal(result[[6]], array(1, c(2,2)))
})
