find_deterministic_relations <- function(M, enforceMinimality=FALSE, exp_index = 1) {
  # Output: the complete list of deterministic relations.
  # The key of the list is the concatenation of the operands with a comma, while the values are the determined variables.
  determined_list <- list()
  
  # Create combinations of candidate columns starting from size 1 to size n.
  for (i in 1:ncol(M)) {
    candidate_operands <- combn(1:ncol(M), i)
    # Enumerate the combinations:
    for (c in 1:ncol(candidate_operands)){
      candidate_operand <- candidate_operands[,c]
      #cat("\nCandidate columns: ", candidate_operand, "\n")
      determined <- array(TRUE, c(1, ncol(M)))
      
      # Unique values of the candidate columns:
      submatrix_det <- M[, candidate_operand]
      vals_operand <- unique(submatrix_det)
      
      # A trick to convert a vector to a matrix.
      if (is.null(nrow(vals_operand))) {
        vals_operand <- matrix(vals_operand, ncol=1)
        submatrix_det <- matrix(submatrix_det, ncol=1)
      } else {
        if (nrow(vals_operand) == nrow(M)) {
          # This combination fully determines the regime and therefore all other variables.
          next
        }
      }
      for (k in 1:nrow(vals_operand)) { 
        # For each value of this combination, get the rows in which it appears.
        val <- vals_operand[k,]
        row_val <- 1:nrow(submatrix_det)
        row_val <- row_val[sapply(1:nrow(submatrix_det), function(x) all(paste(val, collapse = '') %in% paste(submatrix_det[x, ], collapse = '')))]
        val_submatrix <- M[row_val,]
        
        # If there is no ncol, then it's a vector and it's trivially determined.
        if (!is.null(ncol(val_submatrix))){
          # For each column in the rows that match these values:
          for (j in 1:ncol(val_submatrix)){           
            unique_vals <- unique(val_submatrix[, j])
            # If there is more than one value in these cells, the variable is not fully determined.
            if (length(unique_vals) != 1 ) {
              #cat(j, "is not determined by ", candidate_operand, "\n")
              determined[j] <- FALSE
            }
          }
        }
      }
      candidate_operand_label <- paste(candidate_operand, collapse=",")
      # The set of variables minimally determined by this specific candidate operand combination.
      minimally_determined <- c()
      for (k in 1:ncol(M)){
        if (determined[k]) {
          # Check for minimality.
          isMinimal <- TRUE
          if (enforceMinimality && length(candidate_operand) > 1 ) {
            # Check all possible subsets of the candidate_operand combination.
            for (l in 1:(length(candidate_operand)-1)){
              candidate_operand_subsets <- combn(candidate_operand, l)
              for (m in 1:ncol(candidate_operand_subsets)){ 
                candidate_operand_subset <- candidate_operand_subsets[,m]
                # Create the label for the subset.
                candidate_operand_subset_label <- paste(candidate_operand_subset, collapse=",")
                # Get all the variables that are determined by this subset.
                determined_by_subset <- determined_list[[candidate_operand_subset_label]]
                determined_by_subset <- unlist(determined_by_subset)
                if (k %in% determined_by_subset){
#                   if (k != candidate_operand_subset) {
#                     cat("\n", k, ":", candidate_operand, ", determined by subset ",candidate_operand_subset_label, "\n")
#                     print(determined_by_subset)
#                     cat("\n")
#                   }
                  isMinimal <- FALSE
                }
              }
            }
          }
          if (isMinimal) {
            if (k != candidate_operand) {
		cat(k + exp_index - 1, "is determined by ", candidate_operand + exp_index - 1, "\n")
            }
            minimally_determined <- c(minimally_determined, k)
          }
        }
      }
      determined_list[[candidate_operand_label]] <- minimally_determined
    }
  }
  
  determined_list
}


getDet <- function(given=c(4,10,11), determined_list) {
  V <- sort(given)
  V_label <- paste(V, collapse=",")
  determined_list[V_label]
}
