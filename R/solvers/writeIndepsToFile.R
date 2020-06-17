writeIndepsToFile <- function(n, tested_independences,write_function, write_data, 
                              indFilenameFullPath = "../tmp/pipeline.ind", debugOutput=TRUE) {    
  # Write out independences.
  if (debugOutput) {
    debugIndName <- paste(indFilenameFullPath, ".debug.txt", sep="")
    debugInd <- file(debugIndName, "w")
  }
  
  indFile <- file(indFilenameFullPath, "w")
  
  cat('node(1..', n, ').\n', sep='', file = indFile)
  cat('%independences and dependences\n', file = indFile)
  
  for (tested_independence in tested_independences ) {  
    write_function(tested_independence, write_data, indFile)
    if (debugOutput) {
      x <- min(tested_independence$vars)
      y <- max(tested_independence$vars)
      C <- tested_independence$C
      J <- tested_independence$J
      
      if (!is.null(tested_independence$names)){
        x <- tested_independence$names[x]
        y <- tested_independence$names[y]
        C <- tested_independence$names[C]
        J <- tested_independence$names[J]
      }
      
      write_constraints.human_readable(tested_independence$independent, x, y, C, J,
                                         tested_independence$w, tested_independence$p, debugInd) 
    }
  }
  close(indFile)
  
  if (debugOutput) {
    close(debugInd)
  }
  tested_independences
}