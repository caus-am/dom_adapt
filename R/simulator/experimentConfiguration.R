experimentConfiguration <- function(n, type) {
  # Creates the exp configuration, describing the interventions per experiment.
  # columns index variables, rows index experiments.
  # 0 = passive observation, 1 = intervention, -1 = unobserved.
  if ( type == 'single') {
    nexp <- n+1
    # In each experiment there is one intervention on each variable,
    # except the last experiment which has no interventions.
    E <- rbind(diag(n),rep(0,n))
  } else if ( type == 'random' ) {
    nexp<-sample(n,1)
    E <- array(NA, c(nexp,n))
    for ( i in 1:nexp) {
      E[i,] <- dec.to.bin(sample(2^n-1,1), n )
    }
  } else if ( type == 'passive') {
    # Only observations.
    E <- array(0, c(1,n)) 
  } else {
    cat('Bad type of experiments given!\n')
    return(NA)
  }   
  E  
}