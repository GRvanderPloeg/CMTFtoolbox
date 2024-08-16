#' Advanced coupled matrix and tensor factorizations
#'
#' @inherit cmtf_opt
#' @param alpha Scalar penalizing the components to be norm 1 (default 1).
#' @param beta Vector of penalty values for each dataset, penalizing the lambda terms (default 1e-3).
#' @param epsilon Scalar value to make it possible to compute the partial derivatives of lambda (default 1e-8).
#'
#' @return output, see [mize::mize()]
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
acmtf_opt = function(Z, numComponents, initialization="random", alpha=1, beta=rep(1e-3, length(Z$object)), epsilon=1e-8, maxit=2500, tol=1e-6, nstart=1, numCores=1){
  numModes = max(unlist(Z$modes))
  numDatasets = length(Z$object)

  # Create nstart models, in parallel if requested
  if((numCores > 1) & (nstart > 1)){
    cl = parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)
    models = foreach::foreach(i=1:nstart) %dopar% {
      opt=list("fn"=function(x){return(acmtf_fun(x,Z))}, "gr"=function(x){return(acmtf_gradient(x,Z))})
      init = CMTFtoolbox::initializeACMTF(Z, numComponents, initialization, output="vect")
      #model = optimx::Rcgmin(init, function(x){return(acmtf_fun(x,Z))}, function(x){return(acmtf_gradient(x,Z))}, control=control)
      model = mize::mize(init, opt, max_iter=maxit, step_tol=tol, method="CG", line_search="MT")
    }
    parallel::stopCluster(cl)
  } else{
    models = list()
    for(i in 1:nstart){
      opt=list("fn"=function(x){return(acmtf_fun(x,Z))}, "gr"=function(x){return(acmtf_gradient(x,Z))})
      init = initializeACMTF(Z, numComponents, initialization, output="vect")
      #models[[i]] = optimx::Rcgmin(init, function(x){return(acmtf_fun(x,Z))}, function(x){return(acmtf_gradient(x,Z))}, control=control)
      models[[i]] = mize::mize(init, opt, max_iter=maxit, step_tol=tol, method="CG", line_search="MT")
    }
  }

  # Find the best model
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
