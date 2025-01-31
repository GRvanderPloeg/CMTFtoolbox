#' Cross-validation of ACMTF-R by jack-knifing
#'
#' @inherit acmtfr_opt
#' @param maxNumComponents Maximum number of components to investigate (default 5).
#' @param nstart Number of randomly initialized models to create per number of components (default 5).
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
ncrossreg = function(Z, Y, maxNumComponents=5, alpha=1, beta=rep(1e-3, length(Z$object)), epsilon=1e-8, pi=0.5, method="CG", cg_update="HS", line_search="MT", max_iter=10000, max_fn=10000, abs_tol=1e-10, rel_tol=1e-10, grad_tol=1e-10, nstart=5, numCores=1){

  numBlocks = length(Z$object)
  numFolds = Z$sizes[1]

  # Run jack-knifed models to determine RMSE of Y

  # Define settings to run in parallel
  settings = expand.grid(f = 1:maxNumComponents, i = seq_len(numFolds))
  colnames(settings) = c("numComponents", "removeIndex")

  # Pre-allocate memory
  Ypred = rep(NA, nrow(settings))
  Yreal = rep(NA, nrow(settings))

  for(i in 1:nrow(settings)){
    numComponents = settings[i,1]
    removeIndex = settings[i,2]

    Xtrain = lapply(Z$object, function(x){x@data[-removeIndex,,]})
    Ztrain = setupCMTFdata(Xtrain, Z$modes)
    Ytrain = Y[-removeIndex]
    Xtest = lapply(Z$object, function(x){x@data[removeIndex,,]})

    model = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain, numComponents=numComponents, alpha=alpha, beta=beta, epsilon=epsilon, pi=pi, method=method, cg_update=cg_update, line_search=line_search, max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, nstart=nstart, numCores=numCores)
    Ypred[i] = npred(model, Xtest, Ztrain)
    Yreal[i] = Y[removeIndex]
  }

  varExpX = matrix(NA, nrow=maxNumComponents, ncol=numBlocks)
  varExpY = rep(NA, nrow=maxNumComponents)
  for(i in 1:maxNumComponents){
    model = CMTFtoolbox::acmtfr_opt(Z, Y, numComponents=numComponents, alpha=alpha, beta=beta, epsilon=epsilon, pi=pi, method=method, cg_update=cg_update, line_search=line_search, max_iter=max_iter, max_fn=max_fn, abs_tol=abs_tol, rel_tol=rel_tol, grad_tol=grad_tol, nstart=nstart, numCores=numCores)
    varExpX[i,] = model$varExp
    varExpY[i] = model$varExpY
  }

  # Prepare plot data
  varExpData = cbind(1:maxNumComponents, varExpX, varExpY) %>%
    dplyr::as_tibble()
  colnames(varExpData) = c("numComponents", paste0("X", 1:numBlocks), "Y")

  RMSEdata = cbind(settings, Yreal, Ypred) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(RMSE = (Yreal - Ypred)^2 / max(removeIndex)) %>%
    dplyr::group_by(numComponents) %>%
    dplyr::summarise(RMSE = sqrt(sum(RMSE)))

  # Prepare plot output
  a = varExpData %>%
    tidyr::pivot_longer(-numComponents) %>%
    ggplot2::ggplot(ggplot2::aes(x = numComponents,y = value, col = as.factor(name))) +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::labs(x="Number of components", y="Variance explained (%)", color="Block") +
    ggplot2::scale_x_continuous(breaks=varExpData$numComponents)

  b = RMSEdata %>%
    ggplot2::ggplot(ggplot2::aes(x = numComponents, y = RMSE)) +
    ggplot2::geom_path() +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Number of components", y = "RMSE") +
    ggplot2::scale_x_continuous(breaks = RMSEdata$numComponents)

  plot = ggpubr::ggarrange(a,b)
  plots = list(plot, a, b)
  return(list("settings"=settings, "varExp"=varExpData, "RMSE"=RMSEdata, "plots"=plots))
}

# Ugly solution to namespace issues caused by dplyr
RMSE <- NULL
