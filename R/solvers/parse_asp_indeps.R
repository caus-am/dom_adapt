parse_asp_indeps <- function(indPath) {
  indFile <- file(indPath, "r" )
  indLines <- readLines(indFile,-1)
  close(indFile)
  tested_independences <- list()
  intervened_variables <- c()
  
  for (indLine in indLines){
    # Ignore the comments.
    splittedIndLine <- strsplit(indLine,"%")[[1]]
    firstPartLine <- splittedIndLine[1]
    if ( substr(firstPartLine,1,5) == "node(") {
      arguments<-as.numeric(unlist(strsplit(substr(firstPartLine,6,nchar(firstPartLine)-2),'\\.\\.')))
      n <- arguments[2]
      break
    }
  }
  all_vars <- 1:n
  
  for (indLine in indLines){
    # Ignore the comments.
    splittedIndLine <- strsplit(indLine,"%")[[1]]
    firstPartLine <- splittedIndLine[1]
    test_result <- list()
    if ( substr(firstPartLine,1,6) == "indep(") {
      test_result$independent <- TRUE
      arguments<-as.numeric(unlist(strsplit(substr(firstPartLine,7,nchar(firstPartLine)-2),',')))
    } else if (substr(firstPartLine,1,4) == "dep(") {
      test_result$independent <- FALSE
      arguments<-as.numeric(unlist(strsplit(substr(firstPartLine,5,nchar(firstPartLine)-2),',')))
    } else{
      next
    }
    if (arguments[1] < arguments [2]) {
      test_result$vars <-list(arguments[1], arguments[2])
    } else {
      test_result$vars <-list(arguments[2], arguments[1])
    }
    test_result$cset <- arguments[3]
    test_result$jset <- arguments[4]
    test_result$mset <- arguments[5]   
    test_result$C <- all_vars[which(rev(dec.to.bin(test_result$cset, n))==1)]
    
    secondPartsLine <- strsplit(splittedIndLine[2], " ")[[1]]
    for (s in secondPartsLine) {
      if ( substr(s,1,2) == "p=") {
        test_result$p <- as.numeric(substr(s,3,nchar(s)))
      } else if (substr(s,1,2) == "w="){
        test_result$w <- as.numeric(substr(s,3,nchar(s)))
      }
    }
    tested_independences[[length(tested_independences) + 1]] <- test_result
  }
  list(tested_independences=tested_independences, n=n)
}