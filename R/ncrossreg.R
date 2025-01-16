#' Cross-validation of ACMTF-R by jack-knifing
#'
#' @inherit acmtfr_opt
#' @param maxNumComponents Maximum number of components to investigate (default 5).
#' @param nstart Number of randomly initialized models to create per number of components (default 1).
#' @param numCores Number of cores to use (default 1). Setting this to higher than 1 will make the algorithm run in parallel.
#'
#' @return A plot showcasing the variance explained in the X blocks and Y as well as the RMSE of Y.
#' @export
#'
#' @examples
#' set.seed(123)
#' A = array(rnorm(25*2), c(25, 2))
#' B = array(rnorm(100*2), c(100, 2))
#' C = array(rnorm(10*2), c(10, 2))
#' D = array(rnorm(100*2), c(100, 2))
#' E = array(rnorm(10*2), c(10, 2))
#'
#' df1 = reinflateTensor(A, B, C)
#' df2 = reinflateTensor(A, D, E)
#' datasets = list(df1, df2)
#' modes = list(c(1,2,3), c(1,4,5))
#' Z = setupCMTFdata(datasets, modes)
#' Y = matrix(A[,1])
#'
#' # Uses poor settings to make the example work.
#' plot = ncrossreg(Z, Y, maxNumComponents=2, max_iter=2, nstart=2)
ncrossreg = function(Z, Y, maxNumComponents=5, alpha=1, beta=rep(1e-3, length(Z$object)), epsilon=1e-8, pi=0.5, cg_update="HS", line_search="MT", max_iter=10000, max_fn=10000, abs_tol=1e-10, rel_tol=1e-10, grad_tol=1e-10, nstart=5, numCores=1){

  numBlocks = length(Z$object)
  numFolds = Z$sizes[1]

  # Define settings to run in parallel
  settings = expand.grid(f = 1:maxNumComponents, i = seq_len(numFolds))
  settings = cbind(settings, rep(row.names(settings), each=nstart))
  colnames(settings) = c("numComponents", "removeIndex", "replicate")

  # Store models and run for-loop over all settings
  models = list()

  if(numCores > 1){
    cl = parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)

    models = foreach::foreach(i=1:nrow(settings)) %dopar% {
      numComponents = settings[i,1]
      removeIndex = settings[i,2]

      Xtrain = lapply(Z$object, function(x){x@data[-removeIndex,,]})
      Ztrain = setupCMTFdata(Xtrain, Z$modes)
      Ytrain = Y[-removeIndex]
      Xtest = lapply(Z$object, function(x){x@data[removeIndex,,]})

      models[[i]] = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain, numComponents=numComponents, alpha=alpha, beta=beta, epsilon=epsilon, pi=pi, cg_update=cg_update, line_search=line_search, max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, nstart=1)
    }

    parallel::stopCluster(cl)
  } else{
    for(i in 1:nrow(settings)){
      numComponents = settings[i,1]
      removeIndex = settings[i,2]

      Xtrain = lapply(Z$object, function(x){x@data[-removeIndex,,]})
      Ztrain = setupCMTFdata(Xtrain, Z$modes)
      Ytrain = Y[-removeIndex]

      models[[i]] = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain, numComponents=numComponents, alpha=alpha, beta=beta, epsilon=epsilon, pi=pi, cg_update=cg_update, line_search=line_search, max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, nstart=1)
    }
  }

  # Gather model statistics
  varExpX = do.call(rbind, lapply(models, FUN=function(x){x$varExp})) * 100
  varExpY = do.call(rbind, lapply(models, FUN=function(x){x$varExpY})) * 100

  Ypred = rep(NA, nrow(settings))
  Yreal = rep(NA, nrow(settings))
  for(i in 1:nrow(settings)){
    removeIndex = settings[i,2]
    Xtrain = lapply(Z$object, function(x){x@data[-removeIndex,,]})
    Ztrain = setupCMTFdata(Xtrain, Z$modes)
    Xtest = lapply(Z$object, function(x){x@data[removeIndex,,]})
    Ypred[i] = npred(models[[i]], Xtest, Ztrain)
    Yreal[i] = Y[removeIndex]
  }

  # Prepare plot output
  a = cbind(settings, varExpX) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Y = varExpY) %>%
    dplyr::select(-removeIndex,-replicate) %>%
    tidyr::pivot_longer(-c(numComponents)) %>%
    dplyr::group_by(numComponents, name) %>%
    dplyr::summarise(m=mean(value),ymin=min(value), ymax=max(value)) %>%
    ggplot2::ggplot(ggplot2::aes(x=numComponents,y=m,col=as.factor(name))) +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=ymin,ymax=ymax), width=.2) +
    ggplot2::geom_point() +
    ggplot2::labs(x="Number of components", y="Variance explained (%)", color="Block")

  b = cbind(settings, Yreal, Ypred) %>%
    dplyr::as_tibble() %>%
    dplyr::select(-replicate) %>%
    dplyr::mutate(replicate = rep(1:nstart, each=nrow(settings)/nstart)) %>%
    dplyr::mutate(res = Ypred - Yreal) %>%
    dplyr::mutate(res2 = res^2) %>%
    dplyr::mutate(term = res2 / numFolds) %>%
    dplyr::group_by(numComponents, replicate) %>%
    dplyr::summarise(RMSE = sqrt(sum(term))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(numComponents) %>%
    dplyr::summarise(m = mean(RMSE)) %>%
    ggplot2::ggplot(ggplot2::aes(x=as.factor(numComponents),y=m)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(x="Number of components", y="RMSE")

  plot = ggpubr::ggarrange(a,b)
  return(plot)
}

# Ugly solution to namespace issues caused by dplyr
m <- NULL
ymin <- NULL
ymax <- NULL
res <- NULL
res2 <- NULL
term <- NULL
RMSE <- NULL
