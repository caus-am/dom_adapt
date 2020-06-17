learn.asp.short <- function(clingoCmd, clingoInputFiles, n, verbose=0){
  # Initially everything is false:
  C <- array(-1, c(n, n))
  G <- array(-1, c(n, n))
  Ge <- array(-1, c(n, n))
  
  # First put to unknown anything that is true in the union of all models (could be true in some or true in all).
  braveClingo <- paste(clingoCmd, " --opt-mode=optN --enum-mode=brave", clingoInputFiles)
  pipeConnection <- pipe(braveClingo)
  open(pipeConnection)
  while(length(line <- readLines(pipeConnection, n=1)) > 0) {
    for (fact in strsplit(line, " ")[[1]]){
      if (length(grep("UNKNOWN", fact)) > 0) {
        # Timeout
        C <- array(0, c(n, n))
        G <- array(0, c(n, n))
        Ge <- array(0, c(n, n))
        return (list(C=C, G=G, Ge=Ge, Gs=array(0, c(n,n)), objective=-1, timeout=TRUE))
      } else if (length(grep("^causes\\(", fact))>0) {
        vars <- strsplit(gsub("causes\\(", "", gsub("\\)","", fact)), ",")
        C[as.numeric(vars[[1]][2]), as.numeric(vars[[1]][1])] <- 0
      } else if (length(grep("^edge\\(", fact))>0) {
        vars <- strsplit(gsub("edge\\(", "", gsub("\\)","", fact)), ",")
        G[as.numeric(vars[[1]][2]), as.numeric(vars[[1]][1])] <- 0
      } else if (length(grep("^conf\\(", fact))>0) {
        vars <- strsplit(gsub("conf\\(", "", gsub("\\)","", fact)), ",")
        Ge[as.numeric(vars[[1]][2]), as.numeric(vars[[1]][1])] <- 0
        Ge[as.numeric(vars[[1]][1]), as.numeric(vars[[1]][2])] <- 0
      } 
    }
  } #end while reading from stdin.
  close(pipeConnection)
  
  # Second put to true anything that is in the intersection of all models (therefore true in all models).
  cautiousClingo <- paste(clingoCmd, " --opt-mode=optN --enum-mode=cautious", clingoInputFiles)
  pipeConnection <- pipe(cautiousClingo)
  open(pipeConnection)
  while(length(line <- readLines(pipeConnection, n=1)) > 0) {
    for (fact in strsplit(line, " ")[[1]]){
      if (length(grep("UNKNOWN", fact)) > 0) {
        # Timeout
        return (list(C=C, G=G, Ge=Ge, Gs=array(0, c(n,n)), objective=-1, timeout=TRUE))
      } else if (length(grep("^causes\\(", fact))>0) {
        vars <- strsplit(gsub("causes\\(", "", gsub("\\)","", fact)), ",")
        C[as.numeric(vars[[1]][2]), as.numeric(vars[[1]][1])] <- 1
      } else if (length(grep("^edge\\(", fact))>0) {
        vars <- strsplit(gsub("edge\\(", "", gsub("\\)","", fact)), ",")
        G[as.numeric(vars[[1]][2]), as.numeric(vars[[1]][1])] <- 1
      } else if (length(grep("^conf\\(", fact))>0) {
        vars <- strsplit(gsub("conf\\(", "", gsub("\\)","", fact)), ",")
        Ge[as.numeric(vars[[1]][2]), as.numeric(vars[[1]][1])] <- 1
        Ge[as.numeric(vars[[1]][1]), as.numeric(vars[[1]][2])] <- 1
      }
    }
  } #end while reading from stdin.
  close(pipeConnection)
  
  list(C=C, G=G, Ge=Ge, Gs=array(0, c(n,n)), objective=-1, timeout=FALSE)
}