#' Vectorize Fac object
#'
#' @param Fac Fac object from CMTF and ACMTF
#'
#' @return Vectorized Fac object
#' @export
#'
#' @examples
#' set.seed(123)
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*2), c(100, 2))
#' C = array(rnorm(10*2), c(10, 2))
#' D = array(rnorm(100*2), c(100,2))
#' E = array(rnorm(10*2), c(10,2))

#' Fac = list(A, B, C, D, E)
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
#' set.seed(123)
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*2), c(100, 2))
#' C = array(rnorm(10*2), c(10, 2))
#' D = array(rnorm(100*2), c(100,2))
#' E = array(rnorm(10*2), c(10,2))
#'
#' df1 = reinflateTensor(A, B, C)
#' df2 = reinflateTensor(A, D, E)
#' datasets = list(df1, df2)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes, normalize=FALSE)
#'
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

      # 20241119 casting to matrix explicitly to fix Y as a vector corner case
      Fac[[i]] = matrix(Fac[[i]][,sorting], nrow=nrow(Fac[[i]]), ncol=ncol(Fac[[i]]))
    }

    if(ACMTFcase){
      Fac[[numModes+1]] = Fac[[numModes+1]][,sorting]
    }
  }

  # Ensure that matrices are given back, even in the one-component case
  Fac = lapply(Fac, as.matrix)
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
#' set.seed(123)
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*2), c(100, 2))
#' C = array(rnorm(10*2), c(10, 2))
#' D = array(rnorm(100*2), c(100,2))
#' E = array(rnorm(10*2), c(10,2))
#' modes = list(c(1,2,3), c(1,4,5))
#'
#' Fac = list(A, B, C, D, E)
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

#' Runs mize iteratively, allowing for debugging per iteration
#'
#' @inheritParams cmtf_opt
#' @param par Initial parameters
#' @param fg Loss and gradient functions in a list (see mize::mize())
#'
#' @return Model object after convergence
#' @export
#'
#' @examples
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*2), c(100, 2))
#' C = array(rnorm(10*2), c(10, 2))
#' D = array(rnorm(100*2), c(100,2))
#' E = array(rnorm(10*2), c(10,2))
#'
#' df1 = reinflateTensor(A, B, C)
#' df2 = reinflateTensor(A, D, E)
#' datasets = list(df1, df2)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes, normalize=FALSE)
#' numComponents = 1
#' initialization="random"
#'
#' init = initializeCMTF(Z, numComponents, initialization, output="vect")
#' fg = list("fn"=function(x){return(CMTFtoolbox::cmtf_fun(x,Z))},
#'          "gr"=function(x){return(CMTFtoolbox::cmtf_gradient(x,Z))})
#' model = mize_runner(init, fg)
mize_runner = function(par, fg, max_iter=10000, max_fn=10000, abs_tol=1e-8, rel_tol=1e-8, grad_tol=1e-8, cg_update="HS", line_search="MT"){
  init = par
  opt = mize::make_mize(fg=fg, max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, method="CG", cg_update=cg_update, line_search=line_search)
  opt = mize::mize_init(opt=opt, par=par, fg=fg, max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol)
  res = mize::mize_step(opt, par, fg)
  step_info = mize::mize_step_summary(res$opt, res$par, fg=fg)

  all_iterations = init
  all_iterations = cbind(all_iterations, res$par)

  opt_f = fg$fn(par)
  res_f = fg$fn(res$par)
  diff_f = opt_f - res_f

  par = res$par
  opt = res$opt

  while((!mize::check_mize_convergence(step_info)$is_terminated) & (diff_f > abs_tol)){
    res = mize::mize_step(opt, par, fg)

    opt_f = fg$fn(par)
    res_f = fg$fn(res$par)
    diff_f = opt_f - res_f

    step_info = mize::mize_step_summary(res$opt, res$par, fg=fg)
    all_iterations = cbind(all_iterations, res$par)

    par = res$par
    opt = res$opt
    # print(step_info$iter)
    # print(diff_f)
  }

  model = opt
  model$all_iterations = all_iterations
  model$par = par
  model$init = init
  model$f = fg$fn(par)
  model$terminate = step_info$terminate
  return(model)
}
