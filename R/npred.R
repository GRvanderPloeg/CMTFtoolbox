#' Predict Y for new data by projecting the data onto the latent space defined by an ACMTF-R model.
#'
#' @param model ACMTF-R model
#' @param newX List object of new data, where each element corresponds to a block
#' @param Z Original input data used for the model
#' @param sharedMode Shared mode between the blocks (default 1).
#'
#' @return Ypred: the predicted value of Y for the new data
#' @export
#'
#' @examples
#' set.seed(123)
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*2), c(100, 2))
#' C = array(rnorm(10*2), c(10, 2))
#' D = array(rnorm(100*2), c(100, 2))
#' E = array(rnorm(10*2), c(10, 2))
#'
#' df1 = reinflateTensor(A, B, C)
#' df2 = reinflateTensor(A, D, E)
#' datasets = list(df1, df2)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes)
#' Y = matrix(A[,1])

#' # Remove a sample and define
#' i = 1
#' Xtest = lapply(Z$object, function(x){x@data[i,,]})
#' Ytest = Y[i]

#' Xtrain = lapply(Z$object, function(x){x@data[-i,,]})
#' Ytrain = Y[-i]
#' Ztrain = setupCMTFdata(Xtrain, Z$modes)

#' model = acmtfr_opt(Ztrain,Ytrain,2,initialization="random",pi=0, nstart=1, max_iter=10)
#' Ypred = npred(model, Xtest, Ztrain, sharedMode=1)
npred = function(model, newX, Z, sharedMode=1){

  numDatasets = length(model$varExp)
  numComponents = length(model$rho)
  numModes = length(model$Fac) - 1 # the last element contains the lambdas
  Fac = model$Fac

  # Find a projection matrix Z, not to be confused with Z containing the datasets.
  Zproj = list()
  for(p in 1:numDatasets){
    modes = Z$modes[[p]]
    idx = which(modes==sharedMode)
    otherModes = modes[-idx]

    lambdas = Fac[[numModes+1]][p,]
    Zproj[[p]] = multiway::krprod(t(lambdas), multiway::krprod(Fac[[otherModes[2]]], Fac[[otherModes[1]]]))
  }

  # Combine Zproj elements and vectorize newX
  vectZ = do.call(rbind, Zproj)

  if(length(dim(newX[[1]])) > 2){
    vectX = lapply(newX, function(x){t(rTensor::k_unfold(rTensor::as.tensor(x), 1)@data)})
    vectX = do.call(rbind, vectX)
  } else{
    vectX = as.matrix(unlist(lapply(newX, c)))
  }

  # Identify missing values in X, remove those from the calculation
  mask = !is.na(vectX)
  vectZ = vectZ[mask,]
  vectX = vectX[mask,]

  # Project vectX onto the latent space to obtain scores per component
  Zplus = pracma::pinv(vectZ)
  newA = Zplus %*% vectX

  # Predict Y
  if(length(dim(newX[[1]])) > 2){
    Ypred = newA %*% model$rho
  } else{
    Ypred = sum(newA * model$rho)
  }

  return(Ypred)
}
