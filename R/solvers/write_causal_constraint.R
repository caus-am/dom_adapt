write_causal_constraint<-function(trueC, ints, append=TRUE,
                                  output=TRUE, verbose=0,
                                  indFilenameFullPath = "../tmp/pipeline.ind") {
  # Function for writing the causality constraints.
  # trueC true causes matrix
  if (append){
    indFile <- file(indFilenameFullPath, "a")
  } else {
    indFile <- file(indFilenameFullPath, "w")
  }
  
  #cat('%causal relations\n')
  n <- dim(trueC)[1]
  for (i in 1:n) {
    for (j in ints) {
      if (trueC[i,j] == 1) {
        cat(':- not causes(',j,',',i,').\n',sep='', file = indFile)
      } else {
        cat(':- causes(',j,',',i,').\n',sep='', file = indFile)
      }
    }
  }

  close(indFile)
}
