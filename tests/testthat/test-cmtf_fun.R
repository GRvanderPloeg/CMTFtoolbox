test_that("f is nonzero if any solution is found", {
  I = 108
  J = 100
  K = 10
  df = array(rnorm(I*J*K), c(I,J,K))
  datasets = list(df, df)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes)
  result = initializeCMTF(Z, 1, initialization="random")
  f = cmtf_fun(fac_to_vect(result), Z)
  expect_gt(f, 0)
})

test_that("f is zero if the perfect solution is found", {
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*4), c(10, 4))

  df1 = reinflateTensor(A, B[,1:2], C[,1:2])
  df2 = reinflateTensor(A, B[,3:4], C[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4,5))
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)
  fakeResult = list(A, B[,1:2], C[,1:2], B[,3:4], C[,3:4])

  f = cmtf_fun(fac_to_vect(fakeResult), Z)
  expect_equal(f, 0)
})

test_that("an error is thrown for 4-way or more", {
  I = 108
  J = 100
  K = 10
  L = 5
  df = array(rnorm(I*J*K*L), c(I,J,K,L))
  datasets = list(df, df)
  modes = list(c(1,2,3,4), c(1,5,6,7))
  Z = setupCMTFdata(datasets, modes)
  result = initializeCMTF(Z, 1, initialization="random")

  expect_error(cmtf_fun(fac_to_vect(result), Z))
})

test_that("f is zero also in the tensor-matrix case", {
  A = array(rnorm(108*2), c(108, 2))
  B = array(rnorm(100*4), c(100, 4))
  C = array(rnorm(10*2), c(10, 2))

  df1 = reinflateTensor(A, B[,1:2], C)
  df2 = reinflateMatrix(A, B[,3:4])
  datasets = list(df1, df2)
  modes = list(c(1,2,3), c(1,4))
  Z = setupCMTFdata(datasets, modes, normalize=FALSE)
  fakeResult = list(A, B[,1:2], C, B[,3:4])

  f = cmtf_fun(fac_to_vect(fakeResult), Z)
  expect_equal(f, 0)
})
