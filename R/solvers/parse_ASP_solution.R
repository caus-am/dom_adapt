parse_asp_solution<-function(n, solver="clingo", encode=NA,
                         sol_file="tmp/pipeline.ind.clingo",
                         multipleMinima,
                         verbose=1) {
  # Reads a Clingo outputted solutions and makes a model structure out of it.
  # n - number of variables
  # solver - clingo (other solvers were possible at a point)
  # encode - the ASP encoding file used, for some 
  # sol_file - the file in which the solution is
  # Note that the order of the variables in the parsed solution gets inverted when we
  # store them in M$G and M$C (getting "causal ancestor" relationships).
  
  M<-list()
  M$G<-array(NA, c(n,n))
  M$Ge<-array(NA, c(n,n))
  M$objective <- -1 #for the SAT versions
  M$C<-array(NA, c(n,n))
  
  if (!solver %in% learn.asp.available_solvers) 
    stop("Solver ", solver, " not available.")
  
  if (file.exists(sol_file)) {
    solFile <- file( sol_file, "r" )
    sol_lines <- readLines(solFile,-1)
    close(solFile)
    if ( solver == 'clingo') {
      M$sumC<-array(0, c(n,n))
      limit_line <-grep("^TIME LIMIT",sol_lines)
      if (length(limit_line) != 0 ) {  
        return(NA) # returning NA if clingo was stopped with time out
      }
      # Check what kind of results are we reading, the long version or the condensed version.
      long_results <- length(grep("^clingo",sol_lines))
      if (long_results) {
        # Long version of the results.
        M <- parse_long_version(M=M, n=n, sol_lines=sol_lines)
      } else {
        # Short version of the results.
        M$num_solutions <- as.numeric(strsplit(sol_lines[1]," ")[[1]][2])
        M$objective <- as.numeric(strsplit(sol_lines[1]," ")[[1]][1])
        # Parse the causal graph representing the sum of all the models.
        M <- parse_single_solution(M, sol_lines)
        M$sumC <- M$C 
        best_model_matrix <- parse_best_model(n, sol_lines)
        M$bestC <- best_model_matrix
      }
      if (M$num_solutions > 1){
        M$C <- M$sumC - M$num_solutions/2
        
        if (verbose){
          cat("Number of solutions:", M$num_solutions, "\n")
          cat("\nUnnormalized count of causes:\n")
          print(M$sumC)
          cat("\nNormalized count of causes:\n")
          print(M$sumC/M$num_solutions)
          cat("\nMerged causes (1 all agree on an edge, 2 some of them do):\n")
          print(M$C)
        }
      }
    } else {
      M$num_solutions <- 1
      M$C<-array(0, c(n,n))
      M <- parse_single_solution(M, sol_lines)
    }
  } else {
    M$G<-array(0, c(n,n))
    M$Ge<-array(0, c(n,n))
    M$C<-array(0, c(n,n))
  }
  
  M
}

parse_single_solution <- function(P, solution) {
  causes<- P$C
  notcauses<- P$C
  for ( s in solution) { #reading of the graph in the solution
    if ( substr(s,1,5) =="edge(") {
      edgeStrings <- strsplit(s," ")[[1]]
      vars<-as.numeric(unlist(strsplit(substr(edgeStrings[1],6,nchar(edgeStrings[1])-1),',')))
      if (length(edgeStrings) > 1){
        number_of_edges <- as.numeric(edgeStrings[2])
      } else {
        number_of_edges <- 1
      }
      P$G[vars[2],vars[1]]<-number_of_edges
    } else if ( substr(s,1,5) =="conf(") {
      confsStrings <- strsplit(s," ")[[1]]
      vars<-as.numeric(unlist(strsplit(substr(confsStrings[1],6,nchar(confsStrings[1])-1),',')))
      if (length(confsStrings) > 1){
        number_of_confs <- as.numeric(confsStrings[2])
      } else {
        number_of_confs <- 1
      }
      P$Ge[vars[2],vars[1]]<-P$Ge[vars[1],vars[2]]<-number_of_confs
    } else if ( substr(s,1,7) == "causes(") {
      causesStrings <- strsplit(s," ")[[1]]
      causesFact <- gsub("=", "", causesStrings[1])
      vars<-as.numeric(unlist(strsplit(substr(causesFact,8,nchar(causesFact)-1),',')))
      if (length(causesStrings) > 1){
        if (causesStrings[2] =="="){
          number_of_causes <- as.numeric(causesStrings[3])
        } else{
          number_of_causes <- as.numeric(causesStrings[2])
        }
      } else {
        number_of_causes <- 1
      }
      causes[vars[2],vars[1]]<-number_of_causes
    } else if ( substr(s,1,8) == "-causes(" || substr(s,1,10) == "notcauses(") {
      if (substr(s,1,8) == "-causes(") {
        notcausesStrings <- strsplit(s," ")[[1]]
        notcausesFact <- gsub("=", "", notcausesStrings[1])
      } else {
        notcausesStrings <- strsplit(s," = ")[[1]]
        notcausesFact <- notcausesStrings[1]
      }
      if (substr(s,1,10) == "notcauses(") { start <- 11 } else { start <- 9}      
      vars<-as.numeric(unlist(strsplit(substr(notcausesFact,start,nchar(notcausesFact)-1),',')))
      if (length(notcausesStrings) > 1){
        number_of_notcauses <- as.numeric(notcausesStrings[2])
      } else {
        number_of_notcauses <- 1
      }
      notcauses[vars[2],vars[1]]<-number_of_notcauses
    } else if ( substr(s,1,7) == "cindep(") {
      cindepStrings <- strsplit(s," ")[[1]]
      cindepFact <- gsub("=", "", cindepStrings[1])
      vars<-as.numeric(unlist(strsplit(substr(cindepFact,8,nchar(cindepFact)-1),',')))
      if (length(vars) == 3 && length(cindepStrings) > 1){
        number_of_cindeps <- as.numeric(cindepStrings[2])
        if (P$num_solutions != 1  && number_of_cindeps != P$num_solutions) {
          cat("[Warning]: Some models disagree on the cindependences:", s, "\n")
        }
      }
    } else if ( substr(s,1,6) == "indep(") {
      indepStrings <- strsplit(s," ")[[1]]
      indepFact <- gsub("=", "", indepStrings[1])
      vars<-as.numeric(unlist(strsplit(substr(indepFact,7,nchar(indepFact)-1),',')))
      if (length(vars) == 2 && length(indepStrings) > 1){
        number_of_indeps <- as.numeric(indepStrings[2])
        if (P$num_solutions != 1  && number_of_indeps != P$num_solutions) {
          cat("[Warning]: Some models disagree on the independences:", s, "\n")
        }
      }
     }
  }
  P$C <- causes - notcauses
    
  P$G[is.na(P$G)]<-0
  P$Ge[is.na(P$Ge)]<-0
  P
}

parse_best_model <- function(n, solution) {
  # The best model is currently encoded as: 
  # best model:  [0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0], TPR=0.5, FPR=0.5
  best_model_line_index <- grep("^best model:", solution)
  best_model_line <- solution[best_model_line_index]
  best_model_vector_string <- strsplit(strsplit(best_model_line,"\\[")[[1]][2], "\\]")[[1]][1]
  best_model_vector <- as.numeric(strsplit(best_model_vector_string, ",")[[1]])
  best_model_matrix <- matrix(best_model_vector, nrow = n, ncol = n, byrow = FALSE) 
  best_model_matrix
}

parse_long_version <- function(M, n, sol_lines) {
  # Get the compact version of all possible alternatives, with 2 representing an edge
  # that can be 1 in one causal graph and 0 in another.
  solutions <- grep("^Answer",sol_lines)
  M$all_solutions <- list()
  i <- 1
  for (solution in solutions){
    S <- list()
    S$G<-array(0, c(n,n))
    S$Ge<-array(0, c(n,n))
    S$C<-array(0, c(n,n))
    sol <- strsplit(sol_lines[solution + 1],' ')[[1]]
    S <- parse_single_solution(S, sol) 
    # Count duplicates of the current causal graph.
    keepGraph <- TRUE
    for (old_solutions in M$all_solutions) {
      if (length(which(S$C !=  old_solutions$C)) == 0) {
        # The current graph is the same as one of the old ones.
        keepGraph <- FALSE
      }
    }
    if (keepGraph){
      M$sumC <- M$sumC + S$C
      opt_i <- min(grep("^Optimization",sol_lines))
      S$objective<-as.numeric(strsplit(sol_lines[opt_i]," ")[[1]][2])
      M$all_solutions[[i]] <- S
      i <- i + 1
    }
  }
  
  M$G <- M$all_solutions[[1]]$G
  M$Ge <- M$all_solutions[[1]]$Ge
  M$C <- M$all_solutions[[1]]$C
  M$bestC <- M$C
  M$objective <- M$all_solutions[[1]]$objective
  M$num_solutions <- length(M$all_solutions)
  M
}

round_value <- function(softTruth, lowerThresholdFor1 = 0.9, upperThresholdFor0 = 0.1) {
  if ( lowerThresholdFor1 < softTruth  && softTruth < 1) {
    1
  } else if ( 0 < softTruth && softTruth < upperThresholdFor0) {
    0
  } else if (upperThresholdFor0 < softTruth && softTruth < lowerThresholdFor1) {
    2
  } else {
    softTruth
  }
}