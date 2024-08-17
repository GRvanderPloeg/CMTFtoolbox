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
        } else if(length(modes) == 2){
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
