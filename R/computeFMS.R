computeFMS = function(Fac1, Fac2, modes){
  numComponents = ncol(Fac1[[1]])
  numDatasets = length(modes)

  # Recognize which modes are shared between all datasets
  inCommon = intersect(modes[[1]], modes[[2]])
  for(i in 1:numDatasets){
    inCommon = intersect(inCommon, modes[[i]])
  }

  if(length(inCommon) == 0){
    stop("No intersection of modes found.")
  }

  modesToCheck = list()
  for(i in 1:numDatasets){
    modesToCheck[[i]] = setdiff(modes[[i]], inCommon)
  }

  # Calculate FMS per dataset
  FMS_result = rep(0, numDatasets)
  for(i in 1:numDatasets){
    allModesToCheck = modesToCheck[[i]]
    numModes = length(allModesToCheck)
    FMS = 0
    numComparisons = 0

    for(j in 1:numComponents){
      total = 1
      for(k in 1:numModes){
        mode = allModesToCheck[k]
        vect1 = as.matrix(Fac1[[mode]][,j])
        vect2 = as.matrix(Fac2[[mode]][,j])

        #term = abs(t(vect1)*vect2) / (norm(vect1)*norm(vect2))
        term = abs(sum(vect1*vect2)) / (norm(vect1)*norm(vect2))
        total = total * term
      }
      FMS = FMS + total
      numComparisons = numComparisons + 1
    }
  }

  if(numComparisons > 0){
    FMS_result[i] = FMS / numComparisons
  }
  else{
    FMS_result[i] = NA
  }

  return(FMS_result)
}
