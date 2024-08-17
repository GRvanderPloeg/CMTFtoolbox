cmtf_fun = function(x, Z){

  numDatasets = length(Z$object)
  numModes = max(unlist(Z$modes))
  Fac = vect_to_fac(x, Z, sortComponents=FALSE)
  reinflatedBlocks = reinflateFac(Fac, Z, returnAsTensor=TRUE)

  f = 0
  for(p in 1:numDatasets){
    modes = Z$modes[[p]]
    reinflatedBlock = reinflatedBlocks[[p]]
    residuals = Z$object[[p]] - reinflatedBlock
    residuals = Z$missing[[p]] * residuals

    Fnorm = rTensor::fnorm(residuals) # verified to work for matrices too
    f = f + 0.5 * Fnorm^2
  }

  return(f)
}
