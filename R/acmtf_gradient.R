#' Calculate gradient of ACMTF model.
#'
#' @inheritParams acmtf_fun
#'
#' @return Vectorized gradient of the ACMTF model.
#' @export
#'
#' @examples
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*4), c(100, 4))
#' C = array(rnorm(10*4), c(10, 4))
#'
#' df1 = reinflateTensor(A, B[,1:2], C[,1:2])
#' df2 = reinflateTensor(A, B[,3:4], C[,3:4])
#' datasets = list(df1, df2)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes, normalize=FALSE)
#'
#' result = acmtf_opt(Z, 1, max_iter=2) # unoptimized CMTF model
#' f = acmtf_gradient(result$par, Z)
acmtf_gradient = function(x, Z, alpha=1, beta=rep(1e-3, length(Z$object)), epsilon=1e-8){

  numDatasets = length(Z$object)
  numModes = max(unlist(Z$modes))
  Fac = vect_to_fac(x, Z)
  numComponents = ncol(Fac[[1]])
  reinflatedBlocks = reinflateFac(Fac, Z, returnAsTensor=TRUE)
  gradient = list()

  # Gradients per mode stored in a list, will be vectorized at the end.
  for(i in 1:numModes){
    gradient[[i]] = array(0L, dim(Fac[[i]]))

    # Gradient as generated per dataset
    # Note: this is different from CMTF because it multiplies the residuals by the lambdas
    for(p in 1:numDatasets){
      modes = Z$modes[[p]]

      if(i %in% modes){
        idx = which(modes==i)
        otherModes = modes[-idx]

        unfoldedX = rTensor::k_unfold(Z$missing[[p]], idx) * rTensor::k_unfold(Z$object[[p]], idx)
        unfoldedXhat = rTensor::k_unfold(Z$missing[[p]], idx) * rTensor::k_unfold(reinflatedBlocks[[p]], idx)

        if(length(modes) == 3){
          gradientMode = (unfoldedXhat - unfoldedX)@data %*% multiway::krprod(t(Fac[[numModes+1]][p,]), multiway::krprod(Fac[[otherModes[2]]], Fac[[otherModes[1]]]))
        } else if(length(modes) == 2){
          gradientMode = (unfoldedXhat - unfoldedX)@data %*% Fac[[otherModes[1]]] %*% diag(Fac[[numModes+1]][p,])
        }
        else{
          stop(paste0("Number of modes is incorrect for block ", p))
        }

        gradient[[i]] = gradient[[i]] + gradientMode
      }
    }

    # Gradient of norm 1 restriction
    gradient[[i]] = gradient[[i]] + alpha * (Fac[[i]] - removeTwoNormCol(Fac[[i]]))
  }

  # Gradient of the lambdas
  gradient[[numModes+1]] = array(0L, dim(Fac[[numModes+1]]))
  for(i in 1:numDatasets){
    modes = Z$modes[[i]]

    for(j in 1:numComponents){
      lambda_r = Fac[[numModes+1]][i,j]
      residuals = reinflatedBlocks[[i]] - Z$object[[i]]

      if(length(modes) == 3){
        for(k in 1:length(modes)){
          mode = modes[k]
          residuals = rTensor::ttm(residuals, t(as.matrix(Fac[[k]][,j])), k)
        }
        gradient[[numModes+1]][i,j] = residuals@data[1,1,1]
      } else if(length(modes) == 2){
        gradient[[numModes+1]][i,j] = t(Fac[[modes[1]]][,j]) %*% residuals@data %*% Fac[[modes[2]]][,j]
      } else{
        stop(paste0("Number of modes is incorrect for block ", i))
      }

      gradient[[numModes+1]][i,j] = gradient[[numModes+1]][i,j] + ((beta[i]/2) * (lambda_r / (sqrt(lambda_r^2+epsilon))))
    }
  }

  g = fac_to_vect(gradient)
  return(g)
}
