test.classic <- function(vars, C, test_data) {
  # Partial correlation test: 
  #
  # Inputs:
  # vars      - variables used to test the independence
  # C         - conditining set as a list of variable indexes
  # test_data - data for the test:  
  #            test_data$N sample size
  #            test_data$p_threshold p-value threshold used
  #            test_data$Cx covariance matrix of the data
  # Outputs:
  # test_result
  #            test_result$vars sorted variables (standard ass. in the rest of the code)
  #            test_result$C conditioning set
  #            test_result$independent TRUE or FALSE
  #            test_results$p p-value of test
  #            test_results$w weight of the in/dependence
  
  
  test_result <- list()
  test_result$vars <- sort(vars)
  test_result$C <- C
  
  test_result$p <- pcor.indep( test_data$Cx, vars[1], vars[2], given=C, test_data$N );
  test_result$independent <- ( test_result$p > test_data$p_threshold )
  
  # weight is always 1
  test_result$w<-1
  test_result
}

pcor.indep <-function( C, x, y, given=NULL, N ) {
  # Returns a p-value of an independence test.
  #INPUT:
  # C         - Covariance matrix of the data.
  # x,y       - Indexes of variables in question.
  # given     - A vector of variable indexes that are conditioned on.
  # N         - The number of samples.
  #OUTPUT:
  # p-value of the test
  
  cat("Test:", x, " ", y , "|", given, "\n")
  
  
  # Calculating first the partial correlation:
  # If the conditioning set is constant, don't use it.
  if (length(given) > 0 && any(C[c(given), c(given)] == 0)) {
    p <- pcor( c(x,y), C)
  } else {
    p <- pcor( c(x,y,given), C)
  }
  
  #returning a p-value of this test 
  pcor.test(p, length(given), N )$pvalue
}

pcor <- function(u, S) {
  if (det(S[u,u]) != 0) {
    # If anything breaks (computationally singular), try decreasing the tolerance.
    k <- solve(S[u,u], tol=1e-50)
    c <- -k[1,2]/sqrt(k[1,1]*k[2,2])
    if (is.na(c)) { 
      cat("Warning: Determinant is different from zero, but still something is wrong, avoid doing this test:", u, "\n")
    }
  } else {
    cat("Warning: Probably deterministic relation, avoid doing this test:", u, "\n")
    c <- NA
  }
  c
}

pcor.test<-function(r,q,n) {
  df <- n - 2 - q
  tval<- r*sqrt(df)/sqrt(1-r*r)
  pv <- 2* pt(-abs(tval),df)
  list(tval=tval,df=df,pvalue=pv)
}