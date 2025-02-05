#' Cross-validation of ACMTF-R by classical K-fold CV with best-model selection per fold
#'
#' This function runs ACMTF-R with cross-validation. A deterministic K–fold partition
#' is used: the subjects are split in order into `cvFolds` groups. For each fold the
#' training set consists of the other folds and the test set is the current fold.
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
#'                (each CV iteration uses 4/5 for training and 1/5 for testing). Default: 2.
#' @param normalize Normalize the X blocks to frobenius norm 1 (default TRUE).
#' @param normY Normalize Y to a specific value, (default: 1).
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
                      normalize = TRUE,
                      normY = 1,
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
                      cvFolds = 2) {

  ## --- Input Checking and CV Partitioning --- ##
  if (!is.list(Z)) stop("Z must be a list containing 'object', 'sizes', and 'modes'.")
  if (is.null(Z$sizes) || is.null(Z$object) || is.null(Z$modes))
    stop("Z must contain 'sizes', 'object', and 'modes'.")
  numSubjects = Z$sizes[1]
  if (nrow(Y) != numSubjects && length(Y) != numSubjects)
    stop("The number of rows (or length) of Y must equal the number of subjects (Z$sizes[1]).")

  # Create CV folds
  indices = seq_len(numSubjects)
  foldsPartition = split(indices, cut(seq_along(indices), breaks = cvFolds, labels = FALSE))
  uniqueFolds = seq_len(cvFolds)

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

                                      testIdx = foldsPartition[[foldID]]
                                      trainIdx = setdiff(seq_len(numSubjects), testIdx)

                                      ## Prepare X
                                      Xtrain_final = list()
                                      Xtest_final = list()
                                      for(p in 1:length(Z$object)){
                                        Xtrain = rTensor::as.tensor(Z$object[[p]]@data[trainIdx, ,])
                                        Xtest = rTensor::as.tensor(Z$object[[p]]@data[testIdx, ,])

                                        # Centering Xtrain
                                        unfoldedXtrain = rTensor::k_unfold(Xtrain, 1)@data
                                        means = colMeans(unfoldedXtrain, na.rm=TRUE)
                                        unfoldedXtrain_cnt = sweep(unfoldedXtrain, 2, means, FUN="-")
                                        Xtrain_cnt = rTensor::k_fold(unfoldedXtrain_cnt, m=1, modes=Xtrain@modes)

                                        # Use the means to center Xtest as well
                                        unfoldedXtest = rTensor::k_unfold(Xtest, 1)@data
                                        unfoldedXtest_cnt = sweep(unfoldedXtest, 2, means, FUN="-")
                                        Xtest_cnt = rTensor::k_fold(unfoldedXtest_cnt, m=1, modes=Xtest@modes)

                                        # Scaling Xtrain
                                        unfoldedXtrain = rTensor::k_unfold(Xtrain_cnt, 2)@data
                                        stds = apply(unfoldedXtrain, 1, function(x){stats::sd(x, na.rm=TRUE)})
                                        unfoldedXtrain_scl = sweep(unfoldedXtrain, 1, stds, FUN="/")
                                        Xtrain_cnt_scl = rTensor::k_fold(unfoldedXtrain_scl, m=2, modes=Xtrain@modes)

                                        # Use the stds to scale Xtest as well
                                        unfoldedXtest = rTensor::k_unfold(Xtest_cnt, 2)@data
                                        unfoldedXtest_scl = sweep(unfoldedXtest, 1, stds, FUN="/")
                                        Xtest_cnt_scl = rTensor::k_fold(unfoldedXtest_scl, m=2, modes=Xtest@modes)

                                        if(normalize){
                                          norm = rTensor::fnorm(Xtrain_cnt_scl)
                                          Xtrain_final[[p]] = Xtrain_cnt_scl@data / norm
                                          Xtest_final[[p]] = Xtest_cnt_scl@data / norm
                                        } else{
                                          Xtrain_final[[p]] = Xtrain_cnt_scl@data
                                          Xtest_final[[p]] = Xtest_cnt_scl@data
                                        }
                                      }

                                      Ztrain = setupCMTFdata(Xtrain_final, Z$modes, normalize=FALSE) # do not normalize again

                                      ## Prepare Y
                                      Ytrain = Y[trainIdx]
                                      Ymean = mean(Ytrain)
                                      Ytrain_cnt = Ytrain - Ymean

                                      Ynorm = norm(Ytrain_cnt, "2")
                                      Ytrain_normalized = Ytrain_cnt / Ynorm

                                      Ytrain_final = Ytrain_normalized * normY

                                      # Each replicate is fitted with a single initialization and one core.
                                      model_fit = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain_final,
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
                                           Ztrain = Ztrain,
                                           model = model_fit,
                                           Xtest = Xtest_final,
                                           Ymean = Ymean,
                                           Ynorm = Ynorm)
                                    }
    parallel::stopCluster(cl)
  } else {
    resultsList = vector("list", nrow(settings))
    for (i in 1:nrow(settings)) {
      currentRow = settings[i, ]
      currentComp = currentRow$numComponents
      foldID = currentRow$fold
      repID = currentRow$replicate

      testIdx = foldsPartition[[foldID]]
      trainIdx = setdiff(seq_len(numSubjects), testIdx)

      ## Prepare X
      Xtrain_final = list()
      Xtest_final = list()
      for(p in 1:length(Z$object)){
        Xtrain = rTensor::as.tensor(Z$object[[p]]@data[trainIdx, ,])
        Xtest = rTensor::as.tensor(Z$object[[p]]@data[testIdx, ,])

        # Centering Xtrain
        unfoldedXtrain = rTensor::k_unfold(Xtrain, 1)@data
        means = colMeans(unfoldedXtrain, na.rm=TRUE)
        unfoldedXtrain_cnt = sweep(unfoldedXtrain, 2, means, FUN="-")
        Xtrain_cnt = rTensor::k_fold(unfoldedXtrain_cnt, m=1, modes=Xtrain@modes)

        # Use the means to center Xtest as well
        unfoldedXtest = rTensor::k_unfold(Xtest, 1)@data
        unfoldedXtest_cnt = sweep(unfoldedXtest, 2, means, FUN="-")
        Xtest_cnt = rTensor::k_fold(unfoldedXtest_cnt, m=1, modes=Xtest@modes)

        # Scaling Xtrain
        unfoldedXtrain = rTensor::k_unfold(Xtrain_cnt, 2)@data
        stds = apply(unfoldedXtrain, 1, function(x){stats::sd(x, na.rm=TRUE)})
        unfoldedXtrain_scl = sweep(unfoldedXtrain, 1, stds, FUN="/")
        Xtrain_cnt_scl = rTensor::k_fold(unfoldedXtrain_scl, m=2, modes=Xtrain@modes)

        # Use the stds to scale Xtest as well
        unfoldedXtest = rTensor::k_unfold(Xtest_cnt, 2)@data
        unfoldedXtest_scl = sweep(unfoldedXtest, 1, stds, FUN="/")
        Xtest_cnt_scl = rTensor::k_fold(unfoldedXtest_scl, m=2, modes=Xtest@modes)

        if(normalize){
          norm = rTensor::fnorm(Xtrain_cnt_scl)
          Xtrain_final[[p]] = Xtrain_cnt_scl@data / norm
          Xtest_final[[p]] = Xtest_cnt_scl@data / norm
        } else{
          Xtrain_final[[p]] = Xtrain_cnt_scl@data
          Xtest_final[[p]] = Xtest_cnt_scl@data
        }
      }

      Ztrain = setupCMTFdata(Xtrain_final, Z$modes, normalize=FALSE) # do not normalize again

      ## Prepare Y
      Ytrain = Y[trainIdx]
      Ymean = mean(Ytrain)
      Ytrain_cnt = Ytrain - Ymean

      Ynorm = norm(Ytrain_cnt, "2")
      Ytrain_normalized = Ytrain_cnt / Ynorm

      Ytrain_final = Ytrain_normalized * normY

      # Fit model
      model = CMTFtoolbox::acmtfr_opt(Ztrain, Ytrain_final,
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

      resultsList[[i]] = list(numComponents = currentComp,
                              fold = foldID,
                              replicate = repID,
                              testIdx = testIdx,
                              Ztrain = Ztrain,
                              model = model,
                              Xtest = Xtest_final,
                              Ymean = Ymean,
                              Ynorm = Ynorm)
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
    Ytest_comp = rep(NA, numSubjects)
    Ypred_comp = rep(NA, numSubjects)
    foldsToUse = uniqueFolds

    for (fold in foldsToUse) {
      key = paste(comp, fold, sep = "_")
      if (!is.null(bestModels[[key]])) {
        bestEntry = bestModels[[key]]
        pred = npred(bestEntry$model, bestEntry$Xtest, bestEntry$Ztrain)
        pred_norm = pred / normY
        pred_cnt = pred_norm * bestEntry$Ynorm
        pred_original = pred_cnt + bestEntry$Ymean

        Ypred_comp[bestEntry$testIdx] = pred_original
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
