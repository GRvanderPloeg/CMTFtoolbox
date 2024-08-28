#' Factor Match Score (FMS) for model selection
#'
#' @inheritParams acmtf_opt
#' @inheritParams setupCMTFdata
#' @param sharedMode Mode that is shared between all blocks, used to remove fibers for numFolds randomly initialized models.
#' @param model Model type to run, either "acmtf" or "cmtf" (default "acmtf").
#' @param minNumComponents Minimum number of components to check (default 1).
#' @param maxNumComponents Maximum number of components to check (default 3).
#' @param numFolds Number of randomly initialized models to create (default 10).
#' @param jackKnife Jack-knife samples instead of removing multiple samples per fold (default FALSE).
#' @param numCores Number of cores to use (default 1). A number higher than 1 will run the process in parallel.
#'
#' @return List containing "FMS" with the resulting pairwise comparisons of all models per number of components and "plot" with an overview plot.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom foreach %dopar%
#'
#' @examples
#' set.seed(123)
#' A = array(rnorm(108*2), c(108, 2))
#' B = array(rnorm(100*4), c(100, 4))
#' C = array(rnorm(10*4), c(10, 4))
#'
#' df1 = reinflateTensor(A, B[,1:2], C[,1:2])
#' df2 = reinflateTensor(A, B[,3:4], C[,3:4])
#' datasets = list(df1, df2)
#' modes = list(c(1,2,3), c(1,4,5))
#'
#' # specific setting to reduce runtime for CRAN
#' result = investigateFMS(datasets, modes, 1, model="acmtf", numFolds=2, rel_tol=1e-4, abs_tol=1e-4)
investigateFMS = function(datasets, modes, sharedMode, model="acmtf", minNumComponents=1, maxNumComponents=3, numFolds=10, jackKnife=FALSE, numCores=1, max_iter=10000, max_fn=100000, rel_tol=1e-8, abs_tol=1e-8){

  numDatasets = length(datasets)
  numModes = max(unlist(modes))
  components = minNumComponents:maxNumComponents
  FMSresult = list()

  # Fix sizes to only state unique dimensions corresponding to indices in modes
  sizes = rep(0, numModes)
  for(i in 1:numModes){
    for(p in 1:numDatasets){
      if(i %in% modes[[p]]){
        sizes[i] = dim(datasets[[p]])[modes[[p]] == i]
      }
    }
  }
  sharedSize = sizes[sharedMode]

  # A number of fibers of the shared mode is excluded based on numFolds.
  # If numFolds is larger than the size of the shared mode, random jack-knifing will be performed instead.
  # Jack-knifing will always be performed if specified.
  removal = max(1, round(sharedSize / numFolds))

  # Prepare datasets for every fold
  allZ = list()
  for(n in 1:numFolds){
    newDatasets = list()

    # Remove fibers from shared mode
    for(p in 1:numDatasets){
      newSize = dim(datasets[[p]])

      df = rTensor::as.tensor(datasets[[p]])
      df_unfolded = rTensor::k_unfold(df, sharedMode)
      if(jackKnife == TRUE){
        newSize[sharedMode] = newSize[sharedMode] - 1
        dfNew_unfolded = df_unfolded[-n,]
      } else{
        newSize[sharedMode] = newSize[sharedMode] - removal
        dfNew_unfolded = df_unfolded[-sample(1:sharedSize, removal),]
      }
      dfNew = rTensor::k_fold(dfNew_unfolded, sharedMode, newSize)@data
      newDatasets[[p]] = dfNew
    }
    allZ[[n]] = setupCMTFdata(newDatasets, modes, normalize=TRUE)
  }

  # Run models (in parallel if requested)
  for(r in 1:length(components)){
    numComponents = components[r]

    if(model == "acmtf"){
      if((numCores > 1) & (numFolds > 1)){
        cl = parallel::makeCluster(numCores)
        doParallel::registerDoParallel(cl)
        models = foreach::foreach(n=1:numFolds) %dopar% {
          result = acmtf_opt(allZ[[n]], numComponents, initialization="nvec", nstart=1, numCores=1, max_iter=max_iter, max_fn=max_fn, rel_tol=rel_tol, abs_tol=abs_tol)
        }
        parallel::stopCluster(cl)
      } else{
        models = list()
        for(n in 1:numFolds){
          models[[n]] = acmtf_opt(allZ[[n]], numComponents, initialization="nvec", nstart=1, numCores=1, max_iter=max_iter, max_fn=max_fn, rel_tol=rel_tol, abs_tol=abs_tol)
        }
      }
    } else if(model == "cmtf"){
      if((numCores > 1) & (numFolds > 1)){
        cl = parallel::makeCluster(numCores)
        doParallel::registerDoParallel(cl)
        models = foreach::foreach(n=1:numFolds) %dopar% {
          result = cmtf_opt(allZ[[n]], numComponents, initialization="nvec", nstart=1, numCores=1, max_iter=max_iter, max_fn=max_fn, rel_tol=rel_tol, abs_tol=abs_tol)
        }
        parallel::stopCluster(cl)
      } else{
        models = list()
        for(n in 1:numFolds){
          models[[n]] = cmtf_opt(allZ[[n]], numComponents, initialization="nvec", nstart=1, numCores=1, max_iter=max_iter, max_fn=max_fn, rel_tol=rel_tol, abs_tol=abs_tol)
        }
      }
    } else{
      stop("The 'model' argument should be either 'acmtf' or 'cmtf'.")
    }

    FMS = list()
    iterator = 1
    for(i in 1:(numFolds-1)){
      for(j in (i+1):numFolds){
        FMS[[iterator]] = as.data.frame(t(as.matrix(computeFMS(models[[i]]$Fac, models[[j]]$Fac, modes))))
        iterator = iterator + 1
      }
    }
    df = dplyr::bind_rows(FMS) %>% dplyr::as_tibble()
    colnames(df) = paste0("X", 1:numDatasets)
    FMSresult[[r]] = df %>% dplyr::mutate(numComponents = numComponents)
  }

  plottableData = dplyr::bind_rows(FMSresult) %>% tidyr::pivot_longer(-numComponents)
  plot = plottableData %>% dplyr::mutate(numComponents = as.factor(numComponents)) %>% ggplot2::ggplot(ggplot2::aes(x=numComponents,y=value)) +
    ggplot2::facet_wrap(~name) +
    ggplot2::geom_boxplot() +
    ggplot2::xlab("Number of components") +
    ggplot2::ylab("Factor Match Score")

  return(list("FMS"=FMSresult, "plot"=plot))
}

# Ugly solution to namespace issues caused by dplyr
numComponents <- NULL
value <- NULL
name <- NULL
