test_that("setupCMTFdata works with N blocks", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df, df, df, df)
  modes = list(c(1,2,3), c(1,4,5), c(1,6,7), c(1,8,9), c(1,10,11))
  expect_no_error(setupCMTFdata(datasets, modes))
})

test_that("setupCMTFdata throws errors if you don't supply modes", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  expect_error(setupCMTFdata(datasets))
})

test_that("setupCMTFdata identifies the sizes of the datasets correctly", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  expect_equal(Z$sizes, c(108,100,10,100,10))
})

test_that("setupCMTFdata creates a correct missing tensor", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  dfMissing = df
  dfMissing[1,1,1] = NA
  datasets = list(df, dfMissing)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  expect_equal(Z$missing[[2]]@data[1,1,1], 0)
})

test_that("setupCMTFdata sets missing values to zero", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  dfMissing = df
  dfMissing[1,1,1] = NA
  datasets = list(df, dfMissing)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  expect_equal(Z$object[[2]]@data[1,1,1], 0)
})

test_that("setupCMTFdata normalizes correctly when requested", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  expect_equal(rTensor::fnorm(Z$object[[1]]), 1)
})

test_that("setupCMTFdata does not normalize when requested", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)
  expect_equal(rTensor::fnorm(Z$object[[1]]), Z$norms[1])
})

test_that("setupCMTFdata when not normalizing does not change the original data", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)
  expect_equal(df[1,1,1], Z$object[[1]]@data[1,1,1])
})

test_that("fac_to_vect and vect_to_fac work correctly in the CMTF case", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  init = initializeCMTF(Z, 2)

  vect = fac_to_vect(init)
  Fac = vect_to_fac(fac_to_vect(init), Z)
  expect_equal(Fac, init)
})

test_that("fac_to_vect and vect_to_fac work correctly in the ACMTF case", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  init = initializeACMTF(Z, 2)

  vect = fac_to_vect(init)
  Fac = vect_to_fac(fac_to_vect(init), Z)
  expect_equal(Fac, init)
})

test_that("vect_to_fac sorts components correctly", {
  set.seed(123)
  A = array(rnorm(108*5), c(108, 5))
  B = array(rnorm(100*5), c(100, 5))
  C = array(rnorm(10*5), c(10, 5))
  D = array(rnorm(100*5), c(100,5))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateMatrix(A, D)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4))
  Z = setupCMTFdata(datasets, modes)

  result1 = cmtf_opt(Z, 5, initialization="nvec", max_iter=5, sortComponents=FALSE)
  result2 = cmtf_opt(Z, 5, initialization="nvec", max_iter=5, sortComponents=TRUE)
  inputFac = list(A,B,C,D)
})

test_that("removeTwoNormCol indeed removed the two-norm", {
  df = array(rnorm(108,2), c(108,2))
  result = removeTwoNormCol(df)
  expect_equal(norm(result[,1],"2"), 1)
})

test_that("normalizeFac throws no errors", {
  A = array(rnorm(108*2), c(108,2))
  B = array(rnorm(100*2), c(100,2))
  C = array(rnorm(10*2), c(10,2))
  Fac = list(A,B,C)
  expect_no_error(normalizeFac(Fac))
})

test_that("normalizeFac indeed normalized the factors", {
  A = array(rnorm(108*2), c(108,2))
  B = array(rnorm(100*2), c(100,2))
  C = array(rnorm(10*2), c(10,2))
  Fac = list(A,B,C)
  output = normalizeFac(Fac)
  expect_equal(apply(output$Fac[[1]], 2, function(x){norm(as.matrix(x),"F")}), c(1,1))
})

test_that("normalizeFac output dimensions are the same as the input", {
  A = array(rnorm(108*2), c(108,2))
  B = array(rnorm(100*2), c(100,2))
  C = array(rnorm(10*2), c(10,2))
  Fac = list(A,B,C)
  output = normalizeFac(Fac)
  expect_equal(lapply(Fac, dim), lapply(output$Fac, dim))
})

test_that("calculateVarExp does not throw errors", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*2), c(100, 2))
  C = array(rnorm(10*2), c(10, 2))
  D = array(rnorm(100*2), c(100,2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateMatrix(A, D)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4))
  Z = setupCMTFdata(datasets, modes)

  result = cmtf_opt(Z, 2, initialization="nvec", max_iter=5)
  inputFac = list(A,B,C,D)

  expect_no_error(calculateVarExp(result$Fac, Z))
})

test_that("calculateVarExp has a lower bound of zero", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*2), c(100, 2))
  C = array(rnorm(10*2), c(10, 2))
  D = array(rnorm(100*2), c(100,2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateMatrix(A, D)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4))
  Z = setupCMTFdata(datasets, modes)

  result = cmtf_opt(Z, 2, initialization="nvec", max_iter=5)
  inputFac = list(A,B,C,D)

  varExps = calculateVarExp(result$Fac, Z)

  expect_true(all(varExps >= 0))
})

test_that("calculateVarExp has an upper bound of one", {
  set.seed(123)
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*2), c(100, 2))
  C = array(rnorm(10*2), c(10, 2))
  D = array(rnorm(100*2), c(100,2))

  df1 = reinflateTensor(A, B, C)
  df2 = reinflateMatrix(A, D)
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4))
  Z = setupCMTFdata(datasets, modes)

  result = cmtf_opt(Z, 2, initialization="nvec", max_iter=5)
  inputFac = list(A,B,C,D)

  varExps = calculateVarExp(result$Fac, Z)

  expect_true(all(varExps <= 1))
})
