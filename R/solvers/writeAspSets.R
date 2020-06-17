writeAspSets <- function(n, aspSetsFullPath, tested_independences) {
  # Writes down the set definitions for the older encoding.
  # Predicates and "cset", "jset" and "ismember".
  aspSetsFile <- file(aspSetsFullPath, "w")
  
  cat('node(1..', n, ').\n', sep='', file = aspSetsFile)
  # Only write the necessary sets..
  jset_written <- rep(FALSE, 2^n) 
  cset_written <- rep(FALSE, 2^n)
  
  for ( indep in tested_independences ) {
    if ( !cset_written[indep$cset+1] ) {
      # Writing cset only when it is needed
      cat('cset(',indep$cset,'). ',sep='', file = aspSetsFile)
      for ( el in indep$C) {
        cat('ismember(',indep$cset,',',el,'). ',sep='', file = aspSetsFile)        
      }
      cat('\n', file = aspSetsFile)
      cset_written[indep$cset+1]<-TRUE
    }
    
    if ( !jset_written[indep$jset+1] ) {
      cat('jset(',indep$jset,'). ',sep='', file = aspSetsFile)      
      for ( el in indep$J) {
        cat('ismember(',indep$jset,',',el,'). ',sep='', file = aspSetsFile)       
      }
      cat('\n', file = aspSetsFile)
      jset_written[indep$jset+1]<-TRUE      
    }
    
  }#for indep
  cat('\n', file = aspSetsFile)
  close(aspSetsFile)
}


