#' Title
#'
#' @inheritParams acmtfr_opt
#' @param numComponents Vector of integer values containing the number of components F to check.
#' @param pis Vector of pi values to check.
#'
#' @return List object containing "models": all models, "settings": pi and F used to create the models.
#' @export
#'
#' @examples
#' set.seed(123)
#'
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*2), c(100, 2))
#' C = array(rnorm(10*2), c(10, 2))
#' D = array(rnorm(100*2), c(100,2))
#' E = array(rnorm(10*2), c(10,2))
#'
#' df1 = reinflateTensor(A, B, C)
#' df2 = reinflateTensor(A, D, E)
#' Y = A[,1]
#' datasets = list(df1, df2)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes, normalize=FALSE)
#'
#' # specific setting to reduce runtime for CRAN
#' result = cv_degeneracy(Z,
#'                        Y,
#'                        numComponents=1:3,
#'                        pis=c(0.50, 0.75),
#'                        rel_tol=1e-5,
#'                        abs_tol=1e-5,
#'                        nstart=1,
#'                        numCores=1)
cv_degeneracy = function(Z,
                         Y,
                         numComponents=1:5,
                         alpha=1,
                         beta=rep(1e-3, length(Z$object)),
                         epsilon=1e-8,
                         pis=seq(0.05, 1, 0.05),
                         mu=1e-6,
                         initialization="random",
                         method="CG",
                         cg_update="HS",
                         line_search="MT",
                         max_iter=10000,
                         max_fn=10000,
                         abs_tol=1e-10,
                         rel_tol=1e-10,
                         grad_tol=1e-10,
                         nstart=100,
                         numCores=parallel::detectCores()){

  l = list("numComponents"=numComponents,"pi"=pis)
  settings = do.call(expand.grid, l)
  settings = do.call(rbind, replicate(nstart, settings, simplify=FALSE))

  # Run model
  if(numCores == 1){
    models = list()
    for(i in 1:nrow(settings)){
      comp=settings[i,1]
      piValue = settings[i,2]
      models[[i]]=CMTFtoolbox::acmtfr_opt(Z,Y,numComponents=comp,initialization=initialization,beta=beta,pi=piValue,mu=mu,method=method,cg_update=cg_update,line_search=line_search,max_iter=max_iter,max_fn=max_fn,abs_tol=abs_tol,rel_tol=rel_tol,grad_tol=grad_tol,nstart=1)
    }
  } else{
    cl = parallel::makeCluster(parallel::detectCores())
    doParallel::registerDoParallel(cl)
    models = foreach::foreach(i=1:nrow(settings)) %dopar% {
      comp=settings[i,1]
      piValue = settings[i,2]
      model=CMTFtoolbox::acmtfr_opt(Z,Y,numComponents=comp,initialization=initialization,beta=beta,pi=piValue,mu=mu,method=method,cg_update=cg_update,line_search=line_search,max_iter=max_iter,max_fn=max_fn,abs_tol=abs_tol,rel_tol=rel_tol,grad_tol=grad_tol,nstart=1)
    }
    parallel::stopCluster(cl)
  }

  TCCfunc = function(A){
    d = ncol(A)
    M = matrix(0L, nrow=d, ncol=d)

    for(i in 1:d){
      for(j in 1:d){
        M[i,j] = multiway::congru(A[,i],A[,j])
      }
    }
    return(max(abs(M-diag(d))))
  }
  TCCs = unlist(lapply(models, FUN=function(x){TCCfunc(x$Fac[[1]])}))
  plot = cbind(settings,TCCs) %>%
    dplyr::as_tibble() %>%
    ggplot2::ggplot(ggplot2::aes(x=as.factor(numComponents),y=TCCs)) +
    ggplot2::facet_wrap(~pi) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_y_log10() +
    ggplot2::xlab("Number of components") +
    ggplot2::ylab("TCC")

  result = list("models" = models,
                "settings" = settings,
                "plot" = plot)

  return(result)

}
