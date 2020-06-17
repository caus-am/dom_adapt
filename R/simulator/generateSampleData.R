generateSampleData <- function(n, # number of observed variables.
                               exconf, # experiment configuration
                               test, # type of test
                               N, # total number of samples
                               M, # true model of data.
                               verbose = 1){
  # Function that generates the sample data used in pipeline.R.
  #
  # Returns D, a list of structures (size of rows of E, experimental settings), composed by:
  # D$e    - binary vector of size n, where 1 means intervened var
  # D$M    - an intervened model, i.e. a model where edges into int. vars are cut
  # If not oracle, also:
  # D$data - sample data, NA if infinite sample data requested
  # Note, if exconf='passive' as default, no interventions mean only one experimental setting.
  E <- experimentConfiguration(n, exconf)
  nexp <- nrow(E)
  if ( is.infinite(N) ) {
    samples <- rep(Inf,nexp)
  } else {
    # Divide samples evenly among the experiments.
    samples <- divideSamples(N, nexp)
  }
  D <- list()
  for ( i in 1:nexp) {
    D[[i]] <- createData(M,E[i,],samples[i])
  }
  D
}