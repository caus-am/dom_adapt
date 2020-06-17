test.logp.missing <- function(vars, C, test_data) {
  # Partial correlation test that outputs a weight of logp-logalpha
  # - vars: variables used to test the independence
  # - C: conditining set as a list of variable indexes
  # - test_data$X: data for the test:  
  # - test_data$p_threshold: p-value threshold used

  test_result<-list()
  test_result$vars<-sort(vars)
  test_result$C<-C
  
  cov <- test_data$Cx
  N  <- test_data$N
  relevant_cols <- test_data$X [, c(vars[1], vars[2], C)]

  if (any(is.na(relevant_cols)) && length(C) == 0) {
    test_result$p <- NA
  } else {
    if (any(is.na(relevant_cols))) {
      bad_cols <- apply(relevant_cols, 2, function(x) {any(is.na(x)) })      
      original_bad_cols <- c(vars[1], vars[2], C)
      bad_cols <- original_bad_cols[bad_cols]
      
      bad_rows <- apply(relevant_cols, 1, function(x) {any(is.na(x)) })
      original_bad_rows <- 1:nrow(test_data$X)
      bad_rows <- original_bad_rows[bad_rows]
  
      found_all <- TRUE
      # Run partial correlation without the variables with constant values, otherwise they spoil the results.
      removeFromC <- c()

      for (j in 1:length(bad_cols)) {
        bad_col <- bad_cols[j]
        bad_col_values <- test_data$X[, bad_col]
        original_bad_rows <- is.na(bad_col_values)
        
        # Find the indicator variable of the set in which this is all missing.
        found <- FALSE
        for (c in 1:length(C)) {
          if (C[c] %in% bad_cols) {
            next
          }
          if (all(test_data$X[original_bad_rows, C[c]] == 1)) {
            found <- TRUE
            cat( "Found indicator variable in conditioning set:", C[c], " for missing value:", bad_col, "\n")
            removeFromC <- c(removeFromC, C[c])
          }
        }
        
        if (!found) {
          found_all <- FALSE
        }
      }
      
      if (found_all) {
        good_data <- test_data$X[-bad_rows, ]
        cov <- cov(good_data)
        N <- nrow(good_data)
        newC <- setdiff(C, removeFromC)
        test_result$p<-pcor.indep( cov, vars[1], vars[2], given=newC, N );
      } else {
        test_result$p <- NA
      }
    } else {
      # Calls pcor.test that tests with the null hypothesis pcor=0
      test_result$p<-pcor.indep( cov, vars[1], vars[2], given=C, N );
    }
  }

  # If p < alpha reject null hypothesis pcor=0,  else accept it.
  # For Gaussian vars, pcor=0 iff cond independent.
  test_result$independent <- ( test_result$p > test_data$p_threshold )
  
  test_result$w <- abs(log(test_result$p) - log(test_data$p_threshold)) 
  
  if (is.infinite(test_result$w)) {
    test_result$w  <- sign(test_result$w) * 1e6
  }
  test_result
}