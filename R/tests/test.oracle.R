test.oracle<-function(vars,C,test_data) {
  #Function for conduction one independence test.
  #vars      - variables used to test the independence
  #C         - conditining set as a list of variable indexes
  #test_data - data for the test, here test_data$M is the true model.
  
  test_result<-list()
  test_result$vars<-sort(vars) #sorting vars, this is important that no test is considered twice.
  test_result$C<-C
  
  #Here using the oracle implemented.
  #J is already handled in function "learn".
  
  reachable<-directed_reachable(vars[1],vars[2],C,J=c(),test_data$M)
  
  #independence is the negation of reachable
  test_result$independent<- !reachable
  
  #probability of independence is 1 or 0
  test_result$p<-1*test_result$independent
  
  #weight is put to 1 here
  test_result$w<-1
  test_result
}