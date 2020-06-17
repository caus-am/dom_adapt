test.BIC<-function(vars,C,test_data,verbose=FALSE) {
  # BIC test: 
  #
  # Inputs:
  # vars      - variables used to test the independence
  # C         - conditining set as a list of variable indexes
  # test_data - data for the test:  
  #            test_data$N sample size
  #            test_data$p_threshold prior probability
  #            test_data$X data
  #
  # Outputs:
  # test_result
  #            test_result$vars sorted variables (standard ass. in the rest of the code)
  #            test_result$C conditioning set
  #            test_result$independent TRUE or FALSE
  #            test_results$p probability of independence
  #            test_results$w weight of the in/dependence
  
  test_result<-list()
  test_result$vars<-vars
  test_result$C<-C
  
  X<-test_data$X
  x<-vars[1]
  y<-vars[2]
  
  #First determining the ML parameters. 
  #First regression coefficients.
  
  #Model 1 dependence
  fit1<-lm(X[,x] ~ 1+X[,c(y,C)] )
  r1<-residuals(fit1)
  
  #Model 2 independence
  if (length(C) == 0) {
    fit2<-lm(X[,x] ~ 1)
  } else {
    fit2<-lm(X[,x] ~ 1+X[,C])
  }
  r2<-residuals(fit2)
  
  
  #the ML mean nad sd can be obtained as samples
  #calculating here the log-likelihoods of data given ML parameters
  l1<-sum( dnorm(r1,mean=mean(r1),sd=sd(r1),log=TRUE) )
  l2<-sum( dnorm(r2,mean=mean(r2),sd=sd(r2),log=TRUE) )
  
  #the k parameters for the BIC, the dependend model will have one more parameter
  k1<-length(coefficients(fit1))
  k2<-length(coefficients(fit2))
  
  #l1 is higher here of course
  #bic1 <- (-2)*l1 + k1* log(length(x))
  #bic2 <- (-2)*l2 + k2* log(length(x))
  
  #BIC estimates of of the probabilities
  pbic1 <- l1 - k1/2*log(nrow(X)) #*test_data$p_threshold
  pbic2 <- l2 - k2/2*log(nrow(X)) #*test_data$p_threshold

  #this is the prior: prior probabilities of first depedence and then independence
  priors<-c(1-test_data$p_threshold,test_data$p_threshold)

  #  here could multiply with a prior
  #  p<- exp( logp )*priors / sum ( exp(logp)*priors )
  # for numerical accuracy 
  # p<- exp( -max(logp) )exp( logp )*priors / exp( -max(logp) )*sum ( exp(logp)*priors )
  # p<- exp( logp-max(logp)+log(priors) ) / *sum ( exp( logp-max(logp) +log(priors)) )
  
  #add the priors to the previous to get posteriors 
  logp<-c(pbic1,pbic2)+log(priors)
  
  #here calculating the probability
  p <- exp(logp - max(logp))
  p <- p/sum(p)

  #if the penalized likelihood of the first is better than that of the second we conclude ind.
  test_result$independent <- ( logp[2] > logp[1] )

  #Previously there was a possibility for quad weight
  #quad<-c( sum( (p - c(1,0))^2 ),  sum( (p - c(0,1))^2 ) )
  
  if ( test_result$independent) {
    #independence
    test_result$w<-(logp[2]-logp[1])
    #test_result$quad<- quad[1]-quad[2] 
    
  } else {
    #dependence
    test_result$w<-logp[1]-logp[2]    
    #test_result$quad<-quad[2]-quad[1]
    
  }
  
  test_result$p <- p[2] #putting in the probability of independence
  test_result
}