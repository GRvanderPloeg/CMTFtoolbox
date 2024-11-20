#' Advanced coupled matrix and tensor factorizations
#'
#' @inherit cmtf_opt
#' @param alpha Scalar penalizing the components to be norm 1 (default 1).
#' @param beta Vector of penalty values for each dataset, penalizing the lambda terms (default 1e-3).
#' @param epsilon Scalar value to make it possible to compute the partial derivatives of lambda (default 1e-8).
#'
#' @return List object, similar to [mize::mize()] output. Includes a Fac object of the model, which is a list of components per mode. Also includes an init object giving the initialized input vectors.
#' @export
#' @importFrom foreach %dopar%
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
#'
#' # specific setting to reduce runtime for CRAN
#' model = acmtf_opt(Z, 1, rel_tol=1e-5, abs_tol=1e-5)
acmtf_opt = function(Z, numComponents, initialization="random", alpha=1, beta=rep(1e-3, length(Z$object)), epsilon=1e-8, cg_update="HS", line_search="MT", max_iter=10000, max_fn=10000, abs_tol=1e-10, rel_tol=1e-10, grad_tol=1e-10, nstart=1, numCores=1, sortComponents=TRUE, allOutput=FALSE){
  numModes = max(unlist(Z$modes))
  numDatasets = length(Z$object)

  # Prepare initialization outside of the main loop to allow it as output
  inits = list()
  for(i in 1:nstart){
    inits[[i]] = initializeACMTF(Z, numComponents, initialization, output="vect")
  }

  # Create nstart models, in parallel if requested
  if((numCores > 1) & (nstart > 1)){
    cl = parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)
    models = foreach::foreach(i=1:nstart) %dopar% {
      opt = list("fn"=function(x){return(CMTFtoolbox::acmtf_fun(x,Z,alpha,beta,epsilon))}, "gr"=function(x){return(CMTFtoolbox::acmtf_gradient(x,Z,alpha,beta,epsilon))})
      print(dim(Z$object[[1]])) # somehow this line fixes a bug where parallel gives an error about "dims cannot be of length 0"
      model = mize::mize(par=inits[[i]], fg=opt, max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, method="CG", cg_update=cg_update, line_search=line_search)
    }
    parallel::stopCluster(cl)
  } else{
    models = list()
    for(i in 1:nstart){
      models[[i]] = mize::mize(par=inits[[i]], fg=list("fn"=function(x){return(acmtf_fun(x,Z,alpha,beta,epsilon))}, "gr"=function(x){return(acmtf_gradient(x,Z,alpha,beta,epsilon))}), max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, method="CG", cg_update=cg_update, line_search=line_search)
    }
  }

  # Attach extra model info
  for(i in 1:nstart){
    models[[i]]$Fac = vect_to_fac(models[[i]]$par, Z, sortComponents=sortComponents)
    models[[i]]$init = vect_to_fac(inits[[i]], Z, sortComponents=sortComponents)
    models[[i]]$varExp = calculateVarExp(models[[i]]$Fac, Z)
    models[[i]]$varExpPerComponent = calcVarExpPerComponent(models[[i]]$Fac, Z)
  }

  # Return all models if specified, otherwise return only the best model
  if(allOutput == TRUE){
    return(models)
  }
  else{
    bestModel = 0
    bestObjective = Inf
    for(i in 1:nstart){
      model = models[[i]]
      if(model$f <= bestObjective){
        bestModel = model
        bestObjective = model$f
      }
    }
    return(bestModel)
  }
}
