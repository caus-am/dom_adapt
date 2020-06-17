objective <- function(M, tested_independences, asp=TRUE, weight="log", verbose=0) {
  # Calculates the objective function value, defined by the weight-function.
  # M      - model for which the cost is calculated
  # tested_independences - list of the independence relations
  # asp   - if asp is true the each weights is multiplied by 1000 and rounded, otherwise not.
  # weight - "log" or "constant"
  # verbose - if >= 1, the relations with suffered cost will be printed
  
  # For a greedy procedure this could be done much faster...
  
  independentVarsSign <- '_||_'
  dependentVarsSign <- '_N_'
  
  value <- 0
  i <- 0 # index of the underlying loop.
  
  for ( indep in tested_independences) { 
    i <- i+1
    # Calculate the relation indep in the model M.
    independent <- !directed_reachable(indep$vars[1], indep$vars[2], # var1 indep var2
                                       indep$C, # given conditioning set C
                                       indep$J, # given intervening set J
                                       M) # in model M
    if ( independent != indep$independent) { 
      # If the most likely independence in the learned model M is 
      # different from the the independence in the true model (indep), suffer some loss
      if ( !asp && weight == "log") {
        # If not ASP and weight is "log", sum w.
        value <- value + indep$w          
      } else if ( weight == "log") {
        # If ASP and weight is "log", the true weight is divided by 1000.
        indep$w <- round(1000*indep$w) 
        value <- value + indep$w  
      } else if ( weight == "binary" || weight=="bin" ||weight == "none" || weight=='constant' ) {
        # For constant weights, add 1 to the loss.
        value <- value + 1        
      }
      if ( !indep$independent ) {
        data_mark <- dependentVarsSign
      } else {
        data_mark <- independentVarsSign 
      }
      
      if (!independent) {
        learned_mark <- dependentVarsSign
      } else{
        learned_mark <- independentVarsSign
      }
      
      if (verbose) {
        cat(i,':',indep$vars[1],'?',indep$vars[2],'|', paste(indep$C,collapse=','),
            '||',paste(indep$J,collapse=','),
            'data:',data_mark,'learned:',learned_mark,'suffered loss',indep$w,'\n')
      }
    }
  }
  value
}