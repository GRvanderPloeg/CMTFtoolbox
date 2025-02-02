#' Cross-validation of ACMTF-R by jack-knifing
#'
#' @inherit acmtfr_opt
#' @param maxNumComponents Maximum number of components to investigate (default 5).
#' @param nstart Number of randomly initialized models to create per number of components (default 5).
#' @param numCores Number of cores to use (default 1). Setting this to higher than 1 will make the algorithm run in parallel.
#'
#' @return A list object containing: RMSE and varExp data
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
ncrossreg = function(Z, Y, maxNumComponents = 5,
                     alpha = 1,
                     beta = rep(1e-3, length(Z$object)),
                     epsilon = 1e-8,
                     pi = 0.5,
                     method = "CG",
                     cg_update = "HS",
                     line_search = "MT",
                     max_iter = 10000,
                     max_fn = 10000,
                     abs_tol = 1e-10,
                     rel_tol = 1e-10,
                     grad_tol = 1e-10,
                     nstart = 5,
                     numCores = 1) {

  numBlocks = length(Z$object)
  numFolds  = Z$sizes[1]

  #### 1. Jack-knife predictions ####
  # Create grid of settings (each row gives a number of components and a fold index)
  settings = expand.grid(numComponents = 1:maxNumComponents,
                          removeIndex   = seq_len(numFolds))
  # Repeat each setting nstart times (so that each setting is run nstart times)
  settings = settings[rep(1:nrow(settings), each = nstart), ]
  rownames(settings) = NULL  # clean up row names

  if (numCores > 1) {
    # Use parallelization via foreach if more than one core is specified:
    cl = parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)

    results = foreach::foreach(i = 1:nrow(settings),
                                .combine = rbind,
                                .packages = c("CMTFtoolbox")) %dopar% {
                                  # Get the current settings
                                  numComponents = as.numeric(settings$numComponents[i])
                                  removeIndex   = as.numeric(settings$removeIndex[i])

                                  # Prepare jack-knifed training and test sets
                                  Xtrain = lapply(Z$object, function(x) { x@data[-removeIndex,,] })
                                  Ztrain = setupCMTFdata(Xtrain, Z$modes)
                                  Ytrain = Y[-removeIndex]
                                  Xtest  = lapply(Z$object, function(x) { x@data[removeIndex,,] })

                                  # Run ACMTF (with nstart = 1)
                                  model = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain,
                                                                   numComponents = numComponents,
                                                                   alpha         = alpha,
                                                                   beta          = beta,
                                                                   epsilon       = epsilon,
                                                                   pi            = pi,
                                                                   method        = method,
                                                                   cg_update     = cg_update,
                                                                   line_search   = line_search,
                                                                   max_iter      = max_iter,
                                                                   max_fn        = max_fn,
                                                                   abs_tol       = abs_tol,
                                                                   rel_tol       = rel_tol,
                                                                   grad_tol      = grad_tol,
                                                                   nstart        = 1,        # changed here
                                                                   numCores      = numCores)
                                  # Compute prediction for the held-out index:
                                  pred = npred(model, Xtest, Ztrain)
                                  real = Y[removeIndex]

                                  c(Yreal = real, Ypred = pred)
                                }
    parallel::stopCluster(cl)
  } else {
    # Use a simple for-loop when numCores = 1:
    results = matrix(NA, nrow = nrow(settings), ncol = 2)
    colnames(results) = c("Yreal", "Ypred")

    for(i in 1:nrow(settings)){
      numComponents = as.numeric(settings$numComponents[i])
      removeIndex   = as.numeric(settings$removeIndex[i])

      # Prepare jack-knifed training and test sets
      Xtrain = lapply(Z$object, function(x) { x@data[-removeIndex,,] })
      Ztrain = setupCMTFdata(Xtrain, Z$modes)
      Ytrain = Y[-removeIndex]
      Xtest  = lapply(Z$object, function(x) { x@data[removeIndex,,] })

      # Run ACMTF (with nstart = 1)
      model = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain,
                                       numComponents = numComponents,
                                       alpha         = alpha,
                                       beta          = beta,
                                       epsilon       = epsilon,
                                       pi            = pi,
                                       method        = method,
                                       cg_update     = cg_update,
                                       line_search   = line_search,
                                       max_iter      = max_iter,
                                       max_fn        = max_fn,
                                       abs_tol       = abs_tol,
                                       rel_tol       = rel_tol,
                                       grad_tol      = grad_tol,
                                       nstart        = 1,
                                       numCores      = numCores)
      # Compute prediction for the held-out index:
      results[i, ] = c(Yreal = Y[removeIndex], Ypred = npred(model, Xtest, Ztrain))
    }
  }

  # Extract results from the jack-knifing:
  Yreal = results[, "Yreal"]
  Ypred = results[, "Ypred"]

  #### 2. Variance-explained for X and Y ####
  varExpX = matrix(NA, nrow = maxNumComponents, ncol = numBlocks)
  varExpY = rep(NA, maxNumComponents)

  for(i in 1:maxNumComponents){
    # Here i is used as the number of components
    model = CMTFtoolbox::acmtfr_opt(Z, Y,
                                     numComponents = i,    # loop counter
                                     alpha         = alpha,
                                     beta          = beta,
                                     epsilon       = epsilon,
                                     pi            = pi,
                                     method        = method,
                                     cg_update     = cg_update,
                                     line_search   = line_search,
                                     max_iter      = max_iter,
                                     max_fn        = max_fn,
                                     abs_tol       = abs_tol,
                                     rel_tol       = rel_tol,
                                     grad_tol      = grad_tol,
                                     nstart        = nstart,
                                     numCores      = numCores)
    varExpX[i,] = model$varExp
    varExpY[i]  = model$varExpY
  }

  #### 3. Prepare data ####
  # Assign unique column names while binding columns
  colnames(varExpX) = paste0("X", 1:numBlocks)
  varExpData = cbind(numComponents = 1:maxNumComponents,
                     varExpX,
                     Y = varExpY) %>%
    dplyr::as_tibble()

  RMSEdata = cbind(settings, Yreal, Ypred) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(RMSE = (Yreal - Ypred)^2 / max(removeIndex)) %>%
    dplyr::group_by(numComponents) %>%
    dplyr::summarise(RMSE = sqrt(sum(RMSE)))

  return(list("varExp"   = varExpData,
              "RMSE"     = RMSEdata))
}

# Ugly solution to namespace issues caused by dplyr
RMSE <- NULL
