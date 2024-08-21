#' Calculate gradient of CMTF model.
#'
#' @inheritParams cmtf_fun
#'
#' @return Vectorized gradient of the CMTF model.
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
#' result = cmtf_opt(Z, 1, max_iter=2) # unoptimized CMTF model
#' f = cmtf_gradient(result$par, Z)
cmtf_gradient = function(x, Z){

  numDatasets = length(Z$object)
  numModes = max(unlist(Z$modes))
  Fac = vect_to_fac(x, Z, sortComponents=FALSE)
  reinflatedBlocks = reinflateFac(Fac, Z, returnAsTensor=TRUE)
  gradient = list()

  # Gradients per mode stored in a list, will be vectorized at the end.
  for(i in 1:numModes){
    gradient[[i]] = array(0L, dim(Fac[[i]]))

    for(p in 1:numDatasets){
      modes = Z$modes[[p]]

      if(i %in% modes){
        idx = which(modes==i)
        otherModes = modes[-idx]

        unfoldedX = rTensor::k_unfold(Z$missing[[p]], idx) * rTensor::k_unfold(Z$object[[p]], idx)
        unfoldedXhat = rTensor::k_unfold(Z$missing[[p]], idx) * rTensor::k_unfold(reinflatedBlocks[[p]], idx)

        if(length(modes) == 3){
          gradientMode = (unfoldedXhat - unfoldedX)@data %*% multiway::krprod(Fac[[otherModes[2]]], Fac[[otherModes[1]]])
        } else if((length(modes) == 2)){
          gradientMode = (unfoldedXhat - unfoldedX)@data %*% Fac[[otherModes[1]]]
        }
        else{
          stop(paste0("Number of modes is incorrect for block ", p))
        }

        gradient[[i]] = gradient[[i]] + gradientMode
      }
    }
  }

  # CMTF does not penalize components to norm one, nor does it have lambdas.
  # Hence these terms also do not appear in the gradient calculation.

  g = fac_to_vect(gradient)
  return(g)
}
