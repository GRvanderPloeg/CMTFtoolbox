acmtf_fun = function(x, Z, alpha=1, beta=rep(1e-3, length(Z$object)), epsilon=1e-8){

  numDatasets = length(Z$object)
  numModes = max(unlist(Z$modes))
  Fac = vect_to_fac(x, Z)
  numComponents = ncol(Fac[[1]])
  reinflatedBlocks = reinflateFac(Fac, Z, returnAsTensor=TRUE)

  # Penalty for fit on X
  f = 0
  for(p in 1:numDatasets){
    modes = Z$modes[[p]]
    reinflatedBlock = reinflatedBlocks[[p]]
    residuals = Z$object[[p]] - reinflatedBlock
    residuals = Z$missing[[p]] * residuals

    Fnorm = rTensor::fnorm(residuals) # verified to work for matrices too
    f = f + 0.5 * Fnorm^2
  }

  # Penalty to make the solution norm 1
  for(i in 1:numModes){
    for(j in 1:numComponents){
      f = f + alpha * (norm(as.matrix(Fac[[i]][,j]), "2")-1)^2
    }
  }

  # Penalty on the lambdas
  for(i in 1:numComponents){
    for(p in 1:numDatasets){
      f = f + beta[p] * sqrt(Fac[[6]][p,i]^2 + epsilon)
    }
  }

  return(f)
}
