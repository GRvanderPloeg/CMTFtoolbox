#' Set up datasets for (A)CMTF input
#'
#' @param datasets List of arrays of datasets. Multi-way and two-way may be combined.
#' @param modes Numbered modes per dataset in a list. Example element 1: 1 2 3 and element 2: 1 4 for the X tensor and Y matrix case with a shared subject mode.
#' @param normalize Boolean specifying if the datasets should be normalized to Frobenium norm 1.
#'
#' Note: this function puts zeroes in positions with missing values.
#' The indices of missing data are conserved in the output.
#'
#' @return Z, a list with "object" listing the datasets, "sizes" with their size, "norms" with their norms and "missing" stating the missing data.
#' @export
#'
#' @examples
#' I = 108                             # subject mode
#' J = 100                             # feature mode
#' K = 10                              # condition mode
#'
#' df = array(rnorm(I*J*K), c(I,J,K))  # Create the data
#' datasets = list(df, df)             # Two hypothetical data blocks
#' modes = list(c(1,2,3), c(1,4,5))    # Shared subject mode, feature and time modes separate
#' Z = setupCMTFdata(datasets, modes)
setupCMTFdata = function(datasets, modes, normalize=TRUE){

  numDatasets = length(datasets)
  numModes = max(unlist(modes))
  missing = list()
  sizes = list()

  tensors = list()
  for(p in 1:numDatasets){

    # Find and set missing values to zero, conserve indices
    missingMask = is.na(datasets[[p]])
    convertedDataset = datasets[[p]]
    convertedDataset[missingMask] = 0
    missing[[p]] = rTensor::as.tensor(array(as.numeric(!missingMask), dim(datasets[[p]])))
    tensors[[p]] = rTensor::as.tensor(convertedDataset)
  }

  # Calculate norms and normalize if requested
  norms = rep(1,numDatasets)
  for(p in 1:numDatasets){
    norms[p] = rTensor::fnorm(tensors[[p]])

    if(normalize == TRUE){
      tensors[[p]] = tensors[[p]] / rTensor::fnorm(tensors[[p]])
    }
  }

  # Fix sizes to only state unique dimensions corresponding to indices in modes
  sizes = rep(0, numModes)
  for(i in 1:numModes){
    for(p in 1:numDatasets){
      if(i %in% modes[[p]]){
        sizes[i] = dim(datasets[[p]])[modes[[p]] == i]
      }
    }
  }

  Z = list("object"=tensors, "modes"=modes, "sizes"=sizes, "norms"=norms, "missing"=missing)
  return(Z)
}

#' Vectorize Fac object
#'
#' @param Fac Fac object from CMTF and ACMTF
#'
#' @return Vectorized Fac object
#' @export
#'
#' @examples
#' A = array(rnorm(108*2), c(108,2))
#' B = array(rnorm(100*2), c(100,2))
#' C = array(rnorm(10*2), c(10,2))
#' Fac = list(A, B, C)
#' v = fac_to_vect(Fac)
fac_to_vect = function(Fac){
  return(unlist(Fac))
}

#' Convert vectorized output of (a)cmtf to a Fac list object with all loadings per mode.
#'
#' @param vect Vectorized output of (a)cmtf
#' @param Z Original Z input object (see [setupCMTFdata]).
#' @param sortComponents Sort the order of the components by variation explained (default FALSE).
#'
#' @return Fac: list object with all loadings in all components per mode, ordered the same way as Z$modes.
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
#' Z = setupCMTFdata(datasets, modes)
#' result = cmtf_opt(Z, 2, initialization="random", max_iter = 2)
#' Fac = vect_to_fac(result$par, Z)
vect_to_fac = function(vect, Z, sortComponents=FALSE){
  numDatasets = length(Z$object)
  numModes = max(unlist(Z$modes))
  numComponents = length(vect) / sum(Z$sizes)

  Fac = list()
  startIdx = 1
  for(i in 1:numModes){
    Fac[[i]] = array(0L, c(Z$sizes[i], numComponents))

    for(r in 1:numComponents){
      endIdx = startIdx + Z$sizes[i] - 1
      Fac[[i]][,r] = vect[startIdx:endIdx]
      startIdx = endIdx + 1
    }
  }

  # If there are values leftover, you must have an ACMTF model
  ACMTFcase = FALSE
  if(endIdx < length(vect)){
    ACMTFcase = TRUE
  }

  # If you have an ACMTF model, add the remaining values as lambdas
  if(ACMTFcase){
    Fac[[numModes+1]] = array(0L, c(numDatasets, numComponents))
    for(r in 1:numComponents){
      endIdx = startIdx + numDatasets - 1
      Fac[[numModes+1]][,r] = vect[startIdx:endIdx]
      startIdx = endIdx + 1
    }
  }

  if(sortComponents == TRUE){

    # Find variance explained per component
    varExpsPerComp = rep(0, numComponents)
    for(i in 1:numComponents){
      compFac = list()
      for(j in 1:numModes){
        compFac[[j]] = Fac[[j]][,i]
      }

      if(ACMTFcase){
        compFac[[numModes+1]] = as.matrix(Fac[[numModes+1]][,i]) # add lambdas
      }
      varExps = calculateVarExp(compFac, Z)
      varExpsPerComp[i] = mean(varExps)
    }

    sorting = sort(varExpsPerComp, decreasing=TRUE, index.return=TRUE)$ix

    # sort Fac
    for(i in 1:numModes){
      Fac[[i]] = Fac[[i]][,sorting]
    }

    if(ACMTFcase){
      Fac[[numModes+1]] = Fac[[numModes+1]][,sorting]
    }
  }

  return(Fac)
}

#' Create a tensor out of a set of matrices similar to a component model.
#'
#' @param A I x N matrix corresponding to loadings in the first mode for N components.
#' @param B J x N matrix corresponding to loadings in the second mode for N components.
#' @param C K x N matrix corresponding to loadings in the third mode for N components.
#'
#' @return M, an I x J x K tensor.
#' @export
#'
#' @examples
#' A = rnorm(108)
#' B = rnorm(100)
#' C = rnorm(10)
#' M = reinflateTensor(A,B,C)
reinflateTensor = function(A, B, C){

  # Try to cast to matrix if the input is different
  if(!methods::is(A, "matrix")){
    A = as.matrix(A)
  }
  if(!methods::is(B, "matrix")){
    B = as.matrix(B)
  }
  if(!methods::is(C, "matrix")){
    C = as.matrix(C)
  }

  M = array(tcrossprod(A, multiway::krprod(C, B)), c(nrow(A), nrow(B), nrow(C)))
  return(M)
}

#' Create a matrix from a matrix of scores and loadings similar to a component model.
#'
#' @param A I x N matrix corresponding to scores for N components.
#' @param B J x N matrix corresponding to loadings for N components.
#'
#' @return M, an I x J matrix.
#' @export
#'
#' @examples
#' A = rnorm(108)
#' B = rnorm(100)
#' M = reinflateMatrix(A,B)
reinflateMatrix = function(A, B){
  M = tcrossprod(A, B)
  return(M)
}

#' Reinflate all datablocks from a model Fac object.
#'
#' Basically a wrapper function for [reinflateTensor()] and [reinflateMatrix()].
#'
#' @param Fac Fac object output from CMTF and ACMTF
#' @param Z Z object as generated by [setupCMTFdata()].
#' @param returnAsTensor Boolean to return data blocks as rTensor tensor objects (default FALSE)
#'
#' @return List of data blocks
#' @export
#'
#' @examples
#' I = 108
#' J = 100
#' K = 10
#' df = array(rnorm(I*J*K), c(I,J,K))
#' datasets = list(df, df)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes)
#' result = cmtf_opt(Z, 1, max_iter=2)
#' reinflateFac(result$Fac, Z)
reinflateFac = function(Fac, Z, returnAsTensor=FALSE){
  Fac = lapply(Fac, as.matrix) # Cast to matrix for correct indexation in the one-component case.
  numDatasets = length(Z$object)
  numModes = max(unlist(Z$modes))
  numComponents = ncol(Fac[[1]])
  reinflatedFac = list()

  # Check for ACMTF case
  ACMTFcase = FALSE
  if(length(Fac) > numModes){
    ACMTFcase = TRUE
  }

  for(p in 1:numDatasets){
    modes = Z$modes[[p]]
    componentsToSum = list()
    reinflatedBlock = array(0L, dim(Z$object[[p]]))

    for(i in 1:numComponents){
      lambda = ifelse(ACMTFcase, Fac[[numModes+1]][p,i], 1) # Check for ACMTF model lambdas, otherwise lambda=1

      if(length(modes) == 3){
        reinflatedBlock = reinflatedBlock + lambda * reinflateTensor(Fac[[modes[1]]][,i], Fac[[modes[2]]][,i], Fac[[modes[3]]][,i])
      } else if(length(modes) == 2){
        reinflatedBlock = reinflatedBlock + lambda * reinflateMatrix(Fac[[modes[1]]][,i], Fac[[modes[2]]][,i])
      } else{
        stop("Reinflation of blocks of higher modes than 3 is not yet implemented.")
      }
    }

    if(returnAsTensor == TRUE){
      reinflatedFac[[p]] = rTensor::as.tensor(reinflatedBlock)
    }
    else{
      reinflatedFac[[p]] = reinflatedBlock
    }
  }

  return(reinflatedFac)
}

#' Remove two-norms column-wise from a matrix
#'
#' @param df Matrix of loadings
#'
#' @return Matrix of loadings where the column-wise 2-norm is 1.
#' @export
#'
#' @examples
#' A = array(rnorm(108*4), c(108,4))
#' Anorm = removeTwoNormCol(A)
removeTwoNormCol = function(df){
  norms = apply(df, 2, function(x){norm(as.matrix(x), "2")})
  result = sweep(df, 2, norms, FUN="/")
  return(result)
}

#' Normalize all vectors in model output Fac object to norm 1.
#'
#' @param Fac List object with all components per mode per item.
#' @param modes List object with modes per dataset (see also [setupCMTFdata()])
#'
#' @return List object of normalized Fac object, the extracted norms per loading vector per component, and the norms per dataset per component.
#' @export
#'
#' @examples
#' A = array(rnorm(108*2), c(108,2))
#' B = array(rnorm(100*2), c(100,2))
#' C = array(rnorm(10*2), c(10,2))
#' D = array(rnorm(100*2), c(100,2))
#' Fac = list(A,B,C, D)
#' modes = list(c(1,2,3), c(1,4))
#' output = normalizeFac(Fac, modes)
normalizeFac = function(Fac, modes){
  numComponents = ncol(Fac[[1]])
  numModes = length(Fac)
  numDatasets = length(modes)
  normalizedFac = list()

  # Find norms per component in each mode
  extractedNorms = array(0L, c(numModes, numComponents))
  for(i in 1:numModes){
    extractedNorms[i,] = apply(Fac[[i]], 2, function(x){norm(as.matrix(x), "F")})
    normalizedFac[[i]] = sweep(Fac[[i]], 2, extractedNorms[i,], FUN="/")
  }

  # Find norms per component for each dataset
  outputNorms = array(1L, c(numDatasets, numComponents))
  for(i in 1:numComponents){
    for(j in 1:numDatasets){
      relevantModes = modes[[j]]
      for(k in 1:length(relevantModes)){
        mode = relevantModes[k]
        outputNorms[j,i] = outputNorms[j,i] * extractedNorms[mode,i]
      }
    }
  }

  return(list("Fac"=normalizedFac, "normsPerDataset"=outputNorms, "normsPerLoading"=extractedNorms))
}

calculateVarExp = function(Fac, Z){
  Fac = lapply(Fac, as.matrix) # protection from the 1-component case
  numModes = max(unlist(Z$modes))
  numDatasets = length(Z$object)
  reinflatedData = reinflateFac(Fac, Z, returnAsTensor=TRUE)

  varExps = rep(0, numDatasets)
  for(i in 1:numDatasets){
    residuals = Z$object[[i]] - reinflatedData[[i]]
    residualsMissing = Z$missing[[i]] * residuals
    varExps[i] = 1 - ((rTensor::fnorm(residualsMissing)^2) / (rTensor::fnorm(Z$missing[[i]] * Z$object[[i]])^2))
  }

  return(varExps)
}

calcVarExpPerComponent = function(Fac, Z){
  Fac = lapply(Fac, as.matrix) # protection from the 1-component case
  numComponents = ncol(Fac[[1]])
  numModes = max(unlist(Z$modes))
  numDatasets = length(Z$object)

  varExpsPerComp = array(0L, c(numDatasets,numComponents))
  for(i in 1:numComponents){
    compFac = list()
    for(j in 1:numModes){
      compFac[[j]] = Fac[[j]][,i]
    }

    if(length(Fac) > numModes){
      compFac[[numModes+1]] = Fac[[numModes+1]][,i]
    }
    varExpsPerComp[,i] = calculateVarExp(compFac, Z)
  }

  return(varExpsPerComp)
}

findSharedModes = function(modes){
  numDatasets = length(modes)
  sharedModes = modes[[1]]

  for(i in 2:numDatasets){
    sharedModes = intersect(sharedModes, modes[[i]])
  }

  if(length(sharedModes) == 0){
    stop("No intersection of modes found.")
  } else{
    return(sharedModes)
  }

}
