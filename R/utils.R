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
      tensors[[p]] = tensors[[p]]
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

fac_to_vect = function(Fac){
  return(unlist(Fac))
}

vect_to_fac = function(vect, Z){
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

  if(endIdx < length(vect)){ # then you have an ACMTF model
    Fac[[numModes+1]] = array(0L, c(numDatasets, numComponents))
    for(r in 1:numComponents){
      endIdx = startIdx + numDatasets - 1
      Fac[[numModes+1]][,r] = vect[startIdx:endIdx]
      startIdx = endIdx + 1
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

reinflateFac = function(Fac, Z, returnAsTensor=FALSE){
  numDatasets = length(Z$object)
  numModes = max(unlist(Z$modes))
  reinflatedFac = list()

  for(p in 1:numDatasets){
    modes = Z$modes[[p]]

    if(length(modes) == 2){
      reinflatedBlock = reinflateMatrix(Fac[[modes[1]]], Fac[[modes[2]]])
    } else if(length(modes) == 3){
      reinflatedBlock = reinflateTensor(Fac[[modes[1]]], Fac[[modes[2]]], Fac[[modes[[3]]]])
    } else{
      stop("Reinflation of blocks of higher modes than 3 is not yet implemented.")
    }

    if(length(Fac) > numModes){ # then we are dealing with an ACMTF model
      reinflatedBlock = reinflatedBlock * Fac[[numModes+1]][p] # multiply with lambda
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

removeTwoNormCol = function(df){
  result = array(0L, dim(df))

  for(i in 1:ncol(df)){
    result[,i] = df[,i] / norm(as.matrix(df[,i]), "2")
  }

  return(result)
}
