#' Coupled matrix and tensor factorizations
#'
#' @inherit mize::mize
#' @param Z Combined dataset and mode object as produced by [setupCMTFdata()].
#' @param numComponents Number of components
#' @param initialization Initialization, either "random" (default) or "nvec" for numComponents components of the concatenated data using svd.
#' @param nstart Number of models to produce (default 1). If set higher than one, the package will return the best fitted model.
#' @param numCores Number of cores to use (default 1). If set higher than one, the package will attempt to run in parallel.
#'
#' @return List object, similar to [mize::mize()] output. Also includes a Fac object of the model, which is a list of components per mode.
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
#' model = cmtf_opt(Z, 1)
cmtf_opt = function(Z, numComponents, initialization="random", cg_update="HS", line_search="MT", max_iter=10000, max_fn=10000, ls_max_fn=20, abs_tol=1e-8, rel_tol=1e-8, grad_tol=1e-8, nstart=1, numCores=1, allOutput=FALSE){
  numModes = max(unlist(Z$modes))
  numDatasets = length(Z$object)

  # Create nstart models, in parallel if requested
  if((numCores > 1) & (nstart > 1)){
    cl = parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)
    models = foreach::foreach(i=1:nstart) %dopar% {
      opt = list("fn"=function(x){return(cmtf_fun(x,Z))}, "gr"=function(x){return(cmtf_gradient(x,Z))})
      init = CMTFtoolbox::initializeCMTF(Z, numComponents, initialization, output="vect")
      model = mize::mize(init, opt, max_iter=max_iter, max_fn=max_fn, ls_max_fn=ls_max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, method="CG", cg_update=cg_update, line_search=line_search)
    }
    parallel::stopCluster(cl)
  } else{
    models = list()
    for(i in 1:nstart){
      opt = list("fn"=function(x){return(cmtf_fun(x,Z))}, "gr"=function(x){return(cmtf_gradient(x,Z))})
      init = initializeCMTF(Z, numComponents, initialization, output="vect")
      models[[i]] = mize::mize(init, opt, max_iter=max_iter, max_fn=max_fn, ls_max_fn=ls_max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, method="CG", cg_update=cg_update, line_search=line_search)
    }
  }

  # Return all models if specified, otherwise return only the best model
  if(allOutput == TRUE){
    return(models)
  }
  else{
    bestModel = 0
    bestObjective = Inf
    for(i in 1:nstart){
      output = models[[i]]
      if(output$f <= bestObjective){
        bestModel = output
        bestObjective = output$f
      }
    }
    bestModel$Fac = vect_to_fac(bestModel$par, Z)

    return(bestModel)
  }
}
