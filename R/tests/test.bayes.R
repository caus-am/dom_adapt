test.bayes<-function(vars, C, test_data) {
  # Function implementing (with a hack) the Bayesian test from Margaritis2009.
  # vars      - variables used to test the independence
  # C         - conditining set as a list of variable indexes
  # test_data - data for the test:  
  #            test_data$N sample size
  #            test_data$p_threshold prior probability of independence
  #            - deprecated: test_data$X the samples used, now use test_data$df rather 
  #            test_data$df: dataframe of test_data$X
  #            test_data$alpha
  
  # This function uses the 'deal' package code for calculating the local score.
  # The function is currently not particularly fast but it is not the bottleneck.
  # It would be relatively easy to reuse some of the calculated local scores again.
  # For faster operation it is just easier to use the BIC approximation in test.BIC
  # which often gives as good of a performance.

  # Sorry a rather ugly hack of the deal code. No easy to use score for linear gaussian 
  # existed for R at the moment of writing.
  
  # deal needs a data frame, only take the variables and the conditioning set
  if (is.null(test_data$df)){
    df <- data.frame(test_data$X[,c(vars,C)])
  } else {
    df <- test_data$df[,c(vars,C)]
  }
  
  # here using the package deal to for the test
  nw<-network(df)
    
  #ALL TEST IN THE PAPER WERE RUN WITH BELOW 
  #OPTIMAL PRIOR alpha 1.5 and p_threshold 0.1
  #prior<-jointprior.mod(nw,1.5,phiprior="bottcher")
  
  prior<-jointprior.mod(nw,test_data$alpha,phiprior="bottcher")
  
  #first put in as parents only the conditioning set
  if (length(C) == 0 ) {
    nw$nodes[[1]]$parents<-c()
  } else {
    nw$nodes[[1]]$parents<-index(3,ncol(df))
  }
  
  #thise are just some spells from deal code
  node <- nw$nodes[[1]]
  node <- cond.node( node, nw, prior )
  node$condposterior <- node$condprior
  
  node$loglik <- 0
  node <- learnnode( node, nw, df, timetrace = FALSE )
  
  #finally the logp of the independent model
  logpind<-node$loglik
  
  #then essentially add the second variable to the parents
  if (length(C) == 0 ) {
    nw$nodes[[1]]$parents<-2
  } else {
    nw$nodes[[1]]$parents<-index(2,ncol(df))
  }
  node <- nw$nodes[[1]]
  node <- cond.node( node, nw, prior )
  node$condposterior <- node$condprior
  node$loglik <- 0
  node <- learnnode( node, nw, df, timetrace = FALSE )
  #and get the logp of the depedent model
  logpdep<-node$loglik
  
  #then add the priors in the log space
  #p_threshold is the prior prob of indep.
  priors<-c(1-test_data$p_threshold,test_data$p_threshold)

  logp<-c(logpdep,logpind)+log(priors)
  
  #probability vector
  p <- exp(logp - max(logp))
  p <- p/sum(p)
  
  test_result<-list()
  test_result$vars<-vars
  test_result$C<-C
  
  test_result$independent <- ( logp[2] > logp[1] )
  
  if ( test_result$independent) {
    #independence
    test_result$w<-logp[2]-logp[1]
  } else {
    #dependence
    test_result$w<-logp[1]-logp[2]    
  }
  
  test_result$p<-p[2]; #putting in the probability of independence
  
  test_result$prob_dep<-p[1];
  
  test_result
}