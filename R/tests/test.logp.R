test.logp <- function(vars, C, test_data) {
  # Partial correlation test that outputs a weight of logp-logalpha
  # - vars: variables used to test the independence
  # - C: conditining set as a list of variable indexes
  # - test_data: data for the test:  
  # - test_data$N: sample size
  # - test_data$p_threshold: p-value threshold used
  # - test_data$Cx: the covariance matrix used  

  test_result<-list()
  test_result$vars<-sort(vars)
  test_result$C<-C
  
  # Calls pcor.test that tests with the null hypothesis pcor=0
  test_result$p<-pcor.indep( test_data$Cx, vars[1], vars[2], given=C, test_data$N );
  # If p < alpha reject null hypothesis pcor=0,  else accept it.
  # For Gaussian vars, pcor=0 iff cond independent.
  test_result$independent <- ( test_result$p > test_data$p_threshold )
  
  test_result$w <- abs(log(test_result$p) - log(test_data$p_threshold)) 
  
  if (is.infinite(test_result$w)) {
    test_result$w  <- sign(test_result$w) * 1e6
  }
  
  test_result
}