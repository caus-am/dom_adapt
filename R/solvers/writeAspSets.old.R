writeAspSets.old <- function(n, aspSetsFullPath) {
  # Writes down the set definitions as an input file for clingo.
  aspSetsFile <- file(aspSetsFullPath, "a")
  cat('node(1..', n, ').\n', sep='', file = aspSetsFile)

  # Write the sets, which are indexed from 0 to 2^n-1
  for ( i in index(0,2^n-1) ) {
    ibin<-rev(dec.to.bin(i,n))
    for ( j in index(1,length(ibin)) ) {
      if ( ibin[j] == 1 ) {
        cat('ismember(',i,',',j,'). ',sep='', file = aspSetsFile)
      } else {
        # Do nothing.
      }
    }
    cat('\n', file = aspSetsFile)
  }
  close(aspSetsFile)
}