ncrossreg = function(Z, Y, minNumComponents=1, maxNumComponents=10){

  Z_old = Z
  Y_old = Y
  Xblocks_old = lapply(Z$object, function(x){x@data})

  # For now, assume that mode 1 is shared across all blocks
  numFolds = Z$sizes[1]

  # For now, only check 1 component
  f = 1

  # Remove part of the data for testing
  Ypred = rep(NA, numFolds)
  i = 1

  Xblocks_new = lapply(Z$object, function(x){x@data[-i,,]})
  Z_new = setupCMTFdata(Xblocks_new, Z$modes)
  Y_new = Y[-i]
  model = CMTFtoolbox::acmtfr_opt(Z_new, Y_new, f, max_iter = 10)

  # Predict new sample
  newSample = lapply(Z$object, function(x){x@data[i,,]})

  #lambdas = Fac[[numModes+1]][p,] # p is data block number
  # Z = multiway::krprod(model$Fac[[3]], model$Fac[[2]])
  # # multiway::krprod(t(lambdas), multiway::krprod(Fac[[otherModes[2]]], Fac[[otherModes[1]]]))
  # Zplus = pracma::pinv(Z)
  # ti = Zplus %*% c(newSample[[1]])
}
