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

test_that("removeTwoNormCol indeed removed the two-norm", {
  df = array(rnorm(108,2), c(108,2))
  result = removeTwoNormCol(df)
  expect_equal(norm(result[,1],"2"), 1)
})
