test.fake <- function (x, y, S, suffStat) {
  # - suffStat$tested_independences: tested_indeps with a given J assignment.
  # Returns something similar to a p-value (for "bayes" it is prob. of indep) and when below a threshold it is independent.
  list_indeps <- suffStat$tested_independences
  vars <- sort(c(x,y))
  
  for (l in list_indeps){
    if (any( vars != l$vars)){
      next
    }
    if (length(S)== 0 && length(l$C) == 0) return (l$p) 
    if (length(S)== 0) next
    if (length(l$C) == 0) next
    if (length(l$C) != length(S)) next
    if (any(S != l$C)) next
    return (l$p)
  }
  
  stop("Asking for a test that was not performed: ", x, ",", y, ", {", paste(S, collapse=","), "}")
}


test.wrapper <- function (x, y, S, suffStat) {
  vars <- sort(c(x,y))
  test_result <- suffStat$test_function(vars, S, suffStat)
  
  cset <- rep(0, suffStat$n)
  cset[S] <- 1
  cset<-bin.to.dec(rev(cset))
  
  test_result$C <- S
  test_result$jset<- 0
  test_result$cset <- cset
  test_result$M <-setdiff((1:suffStat$n), c(vars,S))
  test_result$mset <- getm( vars, S, n=suffStat$n)
  
  indFile <- file(suffStat$indPath, "a")
  write_constraint(test_result, list(weight="log"), indFile)
  close(indFile)
  test_result$p
}