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
#' I = 108
#' J = 100
#' K = 10
#' df = array(rnorm(I*J*K), c(I,J,K))
#' datasets = list(df, df)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes)
#' model = acmtf_opt(Z, 1)
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
      model = mize::mize(par=inits[[i]], fg=list("fn"=function(x){return(acmtf_fun(x,Z))}, "gr"=function(x){return(acmtf_gradient(x,Z))}), max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, method="CG", cg_update=cg_update, line_search=line_search)
    }
    parallel::stopCluster(cl)
  } else{
    models = list()
    for(i in 1:nstart){
      models[[i]] = mize::mize(par=inits[[i]], fg=list("fn"=function(x){return(acmtf_fun(x,Z))}, "gr"=function(x){return(acmtf_gradient(x,Z))}), max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, method="CG", cg_update=cg_update, line_search=line_search)
    }
  }

  # Return all models if specified, otherwise return only the best model
  if(allOutput == TRUE){
    output = list()
    for(i in 1:nstart){
      output[[i]] = models[[i]]
      output[[i]]$Fac = vect_to_fac(models[[i]]$par, Z, sortComponents=sortComponents)
      output[[i]]$init = vect_to_fac(inits[[i]], Z, sortComponents=sortComponents)
    }
    return(output)
  }
  else{
    bestModel = 0
    bestObjective = Inf
    bestInit = NA
    for(i in 1:nstart){
      output = models[[i]]
      if(output$f <= bestObjective){
        bestModel = output
        bestObjective = output$f
        bestInit = inits[[i]]
      }
    }
    bestModel$Fac = vect_to_fac(bestModel$par, Z, sortComponents=sortComponents)
    bestModel$init = vect_to_fac(bestInit, Z, sortComponents=sortComponents)

    return(bestModel)
  }
}
