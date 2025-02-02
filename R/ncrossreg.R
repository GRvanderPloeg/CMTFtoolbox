#' Cross-validation of ACMTF-R by classical K-fold CV (or jackknife) with best-model selection per fold
#'
#' This function runs ACMTF-R with cross-validation. If the number of folds
#' (`cvFolds`) is provided and is less than the number of subjects (i.e. the first
#' dimension given in `Z$sizes[1]`), then a classical, deterministic K–fold partition
#' is used: the subjects are split in order into `cvFolds` groups. For each fold the
#' training set consists of the other folds and the test set is the current fold.
#' If `cvFolds` is `NULL` (the default) or is greater than or equal to the number of subjects,
#' then leave–one–out (jackknife) CV is performed.
#'
#' For each fold and for each number of components, *nstart* models are fitted—each
#' with a single initialization (i.e. using `nstart = 1` and `numCores = 1` in `acmtfr_opt`).
#' The large parallel loop iterates over all combinations of (numComponents, fold, replicate).
#' After the loop, for each (numComponents, fold) the best replicate is selected (using `model$f`)
#' and used to predict on the test set. The predictions for all folds are then combined into
#' a unified prediction vector for RMSE calculation.
#'
#' @inheritParams acmtfr_opt
#' @param maxNumComponents Maximum number of components to investigate (default 5).
#' @param nstart Number of replicate initializations per CV fold (default 5).
#' @param numCores Number of cores to use for the replicate fits (default 1). Each replicate model is
#'                 fit with `numCores = 1` (to avoid nested parallelism) but the replicates are run in parallel.
#' @param cvFolds Number of folds to use in the cross-validation. For example, if `cvFolds`
#'                is 5, then the subjects are deterministically partitioned into 5 groups
#'                (each CV iteration uses 4/5 for training and 1/5 for testing).
#'                If `cvFolds` is `NULL` (the default) or is greater than or equal to the number of subjects,
#'                then leave–one–out (jackknife) CV is performed.
#' @param center Center the data per fold Does not apply to Y.
#' @param centerY Center Y per fold
#' @param scale Scale the data per partition. Does not apply to Y.
#' @param normalize Block scale the datasets to norm 1 per fold.
#'
#' @return A list with two elements:
#'         - **varExp**: a tibble with the variance–explained (for X and Y) per number of components.
#'         - **RMSE**: a tibble with the RMSE (computed over the unified CV prediction vector) per number of components.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' A = array(rnorm(25 * 2), c(25, 2))
#' B = array(rnorm(100 * 2), c(100, 2))
#' C = array(rnorm(10 * 2), c(10, 2))
#' D = array(rnorm(100 * 2), c(100, 2))
#' E = array(rnorm(10 * 2), c(10, 2))
#'
#' df1 = reinflateTensor(A, B, C)
#' df2 = reinflateTensor(A, D, E)
#' datasets = list(df1, df2)
#' modes = list(c(1, 2, 3), c(1, 4, 5))
#' Z = setupCMTFdata(datasets, modes)
#' Y = matrix(A[, 1])
#'
#' # For classical 5-fold CV (deterministic splits) with best-model selection:
#' result = ncrossreg(Z, Y, maxNumComponents = 2, max_iter = 2, nstart = 2, cvFolds = 5)
ncrossreg = function(Z, Y,
                      maxNumComponents = 5,
                      alpha = 1,
                      beta = rep(1e-3, length(Z$object)),
                      epsilon = 1e-8,
                      pi = 0.5,
                      center = TRUE,
                      centerY = TRUE,
                      scale = TRUE,
                      normalize = TRUE,
                      method = "CG",
                      cg_update = "HS",
                      line_search = "MT",
                      max_iter = 10000,
                      max_fn = 10000,
                      abs_tol = 1e-10,
                      rel_tol = 1e-10,
                      grad_tol = 1e-10,
                      nstart = 5,
                      numCores = 1,
                      cvFolds = NULL) {

  ## --- Input Checking and CV Partitioning --- ##
  if (!is.list(Z)) stop("Z must be a list containing 'object', 'sizes', and 'modes'.")
  if (is.null(Z$sizes) || is.null(Z$object) || is.null(Z$modes))
    stop("Z must contain 'sizes', 'object', and 'modes'.")
  numSubjects = Z$sizes[1]
  if (nrow(Y) != numSubjects && length(Y) != numSubjects)
    stop("The number of rows (or length) of Y must equal the number of subjects (Z$sizes[1]).")
  if (is.null(cvFolds)) cvFolds = numSubjects  # default to jackknife if not specified
  if(centerY) Y  = Y - mean(Y)

  if (cvFolds < numSubjects) {
    cvMode = "kfold"
    indices = seq_len(numSubjects)
    foldsPartition = split(indices,
                            cut(seq_along(indices), breaks = cvFolds, labels = FALSE))
    uniqueFolds = seq_len(cvFolds)
  } else {
    cvMode = "jackknife"
    uniqueFolds = seq_len(numSubjects)
    foldsPartition = NULL
  }

  ## --- Create Settings Data Frame ---
  settings = expand.grid(numComponents = 1:maxNumComponents,
                          fold = uniqueFolds,
                          replicate = 1:nstart)
  settings = settings[order(settings$numComponents, settings$fold, settings$replicate), ]

  ## --- Run the Parallel Loop Over All Settings ---
  if (numCores > 1) {
    cl = parallel::makeCluster(numCores)
    doParallel::registerDoParallel(cl)
    resultsList = foreach::foreach(i = 1:nrow(settings),
                                    .packages = c("CMTFtoolbox")) %dopar% {
                                      currentRow = settings[i, ]
                                      currentComp = currentRow$numComponents
                                      foldID = currentRow$fold
                                      repID = currentRow$replicate

                                      if (cvMode == "kfold") {
                                        testIdx = foldsPartition[[foldID]]
                                      } else {
                                        testIdx = foldID
                                      }
                                      trainIdx = setdiff(seq_len(numSubjects), testIdx)

                                      Xtrain = lapply(Z$object, function(x) x@data[trainIdx, ,])
                                      Ztrain = setupCMTFdata(Xtrain, Z$modes, normalize=normalize)
                                      Ytrain = Y[trainIdx]

                                      if(center) Xtrain = lapply(Xtrain, parafac4microbiome::multiwayCenter)
                                      if(scale)  Xtrain = lapply(Xtrain, parafac4microbiome::multiwayScale)
                                      if(centerY) Ytrain = Ytrain - mean(Ytrain)

                                      # Each replicate is fitted with a single initialization and one core.
                                      model_fit = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain,
                                                                           numComponents = currentComp,
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
                                                                           numCores      = 1)
                                      list(numComponents = currentComp,
                                           fold = foldID,
                                           replicate = repID,
                                           testIdx = testIdx,
                                           model = model_fit,
                                           Ztrain = Ztrain)
                                    }
    parallel::stopCluster(cl)
  } else {
    resultsList = vector("list", nrow(settings))
    for (i in 1:nrow(settings)) {
      currentRow = settings[i, ]
      currentComp = currentRow$numComponents
      foldID = currentRow$fold
      repID = currentRow$replicate

      if (cvMode == "kfold") {
        testIdx = foldsPartition[[foldID]]
      } else {
        testIdx = foldID
      }
      trainIdx = setdiff(seq_len(numSubjects), testIdx)

      Xtrain = lapply(Z$object, function(x) x@data[trainIdx, ,])
      Ztrain = setupCMTFdata(Xtrain, Z$modes, normalize=normalize)
      Ytrain = Y[trainIdx]

      if(center) Xtrain = lapply(Xtrain, parafac4microbiome::multiwayCenter)
      if(scale)  Xtrain = lapply(Xtrain, parafac4microbiome::multiwayScale)
      if(centerY) Ytrain = Ytrain - mean(Ytrain)

      resultsList[[i]] = list(numComponents = currentComp,
                               fold = foldID,
                               replicate = repID,
                               testIdx = testIdx,
                               Ztrain = Ztrain,
                               model = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain,
                                                               numComponents = currentComp,
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
                                                               numCores      = 1),
                               Ztrain = Ztrain)
    }
  }

  ## --- Group the Results by (numComponents, fold) and Select the Best Model --- ##
  # Create a grouping key for each result.
  keys = sapply(resultsList, function(x) paste(x$numComponents, x$fold, sep = "_"))
  resultsByGroup = split(resultsList, keys)

  bestModels = lapply(resultsByGroup, function(group) {
    f_vals = sapply(group, function(x) x$model$f)
    group[[which.min(f_vals)]]
  })

  ## --- Assemble Predictions and Compute RMSE for Each Number of Components --- ##
  RMSE_list = rep(NA, maxNumComponents)
  predictionsByComp = vector("list", maxNumComponents)

  for (comp in 1:maxNumComponents) {
    Ypred_comp = rep(NA, numSubjects)
    if (cvMode == "kfold") {
      foldsToUse = uniqueFolds
    } else {
      foldsToUse = seq_len(numSubjects)
    }
    for (fold in foldsToUse) {
      key = paste(comp, fold, sep = "_")
      if (!is.null(bestModels[[key]])) {
        bestEntry = bestModels[[key]]
        testIdx = bestEntry$testIdx
        Xtest = lapply(Z$object, function(x) x@data[testIdx, ,])
        if(center) Xtest  = lapply(Xtest, parafac4microbiome::multiwayCenter)
        if(scale)  Xtest  = lapply(Xtest, parafac4microbiome::multiwayScale)
        pred = npred(bestEntry$model, Xtest, bestEntry$Ztrain)
        Ypred_comp[testIdx] = pred
      }
    }
    predictionsByComp[[comp]] = Ypred_comp
    RMSE_list[comp] = sqrt(mean((Y - Ypred_comp)^2))
  }

  RMSEdata = dplyr::tibble(numComponents = 1:maxNumComponents,
                            RMSE = RMSE_list)

  ## --- Fit Full-Data Models to Compute Variance-Explained --- ##
  varExpX = matrix(NA, nrow = maxNumComponents, ncol = length(Z$object))
  varExpY = rep(NA, maxNumComponents)
  for (i in 1:maxNumComponents) {
    full_model = CMTFtoolbox::acmtfr_opt(Z, Y,
                                          numComponents = i,
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
    varExpX[i, ] = full_model$varExp
    varExpY[i]   = full_model$varExpY
  }
  colnames(varExpX) = paste0("X", 1:length(Z$object))
  varExpData = dplyr::as_tibble(cbind(numComponents = 1:maxNumComponents,
                                       varExpX,
                                       Y = varExpY))

  return(list("varExp" = varExpData,
              "RMSE"   = RMSEdata))
}

# Ugly solution to namespace issues caused by dplyr
RMSE = NULL
