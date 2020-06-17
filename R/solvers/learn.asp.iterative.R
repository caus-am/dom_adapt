learn.asp.iterative <- function(clingoCmd, clingoInputFiles, n, graph_filename_template,
                       tested_independences,
                       currentDir = "./../tmp/",
                       debugOutput=FALSE,
                       parallelise=TRUE,
                       evalDirectCauses=FALSE,
                       returnNothingIfTimeout=FALSE,
                       verbose=1){
  multipleMinimaTempDir <- paste(currentDir, "/multipleMinimaTemp/", sep ="")
  system(paste("mkdir -p", multipleMinimaTempDir))
  
  C <- array(0, c(n,n))
  G <- array(0, c(n,n))
  
  baselineClingoCmd <- paste(clingoCmd, clingoInputFiles)
  if (debugOutput) {
    optCmd <- paste(baselineClingoCmd, "|tee",  paste(multipleMinimaTempDir, "/", graph_filename_template, "_baseline.txt", sep=""))
  } else {
    optCmd <- paste(baselineClingoCmd)
  }
  optResult <- system(optCmd, intern = TRUE)
  
  if (length(grep("INTERRUPTED", optResult))> 0 || length(grep("TIME LIMIT", optResult))> 0){
    if (returnNothingIfTimeout) {
      return (list(C=C, G=array(0, c(n,n)), Ge=array(0, c(n,n)), Gs=array(0, c(n,n)), objective=-1, timeout=TRUE))
    }
    baselineTimeout = TRUE
  } else {
    baselineTimeout = FALSE
  }
  
  if (length(grep("UNSATISFIABLE", optResult))> 0) {
    # The optimization is unsatisfiable.
    return (list(C=array(-Inf, c(n,n)), G=array(0, c(n,n)), Ge=array(0, c(n,n)), Gs=array(0, c(n,n)), objective=-1, timeout=FALSE))
  } else if (length(grep("OPTIMUM FOUND", optResult))> 0){
    optimumLine <- optResult[grep('Optimization : ', optResult)]
    optimum_baseline <- as.numeric(gsub('Optimization : ','', optimumLine))
  } else if (length(grep("SATISFIABLE", optResult))> 0){
    optimumLine <- optResult[grep('Optimization : ', optResult)]
    optimum_baseline <- as.numeric(gsub('Optimization : ','', optimumLine))
  } else if (baselineTimeout) {
    optimum <- 0
  } else {
    stop("Strange behaviour of solver:", optResult)
  }
  
  if (parallelise) {
    scoresC <- foreach (i = 1:length(C)) %dopar% {
      learn.asp.iterative.loop (n, graph_filename_template,  (i-1)%% n + 1 , (i-1)%/% n + 1, predicateToTest= "causes",
                                multipleMinimaTempDir, baselineClingoCmd, optimum_baseline=optimum_baseline, debugOutput=debugOutput, verbose=verbose)
    }
    if (evalDirectCauses) {
      ancestralResults <-  paste(currentDir, "/", graph_filename_template, "_ancestral.txt", sep="")
      sink(ancestralResults)
      for (i in 1:length(scoresC)){
        arg1 <- (i-1)%% n + 1 
        arg2 <- (i-1)%/% n + 1
        if(scoresC[[i]]$score > 0) {
          cat(":- not causes(",arg1,",", arg2,").\n", sep="")
        } else if(scoresC[[i]]$score < 0) {
          cat(":- causes(",arg1,",", arg2,").\n", sep="")
        }   
      }
      sink()
      if (grepl("tom", clingoInputFiles) > 0){
        clingoInputFiles <- gsub("ASP/tom_acyclic_complete.pl", "ASP/new_wmaxsat_acyclic_causes.pl", clingoInputFiles)
        clingoInputFiles <- gsub("ASP/tom_acyclic_maxschedule_nodep.pl", "ASP/new_wmaxsat_acyclic_causes.pl", clingoInputFiles)
        aspSetsFullPathDirect <-  paste(currentDir, "/",graph_filename_template, "_direct.pre.asp", sep="")
        writeAspSets(n=n, aspSetsFullPath = aspSetsFullPathDirect, tested_independences=tested_independences)
        build_tree(n=n, aspSetsFullPath = aspSetsFullPathDirect, tested_independences=tested_independences)   
        baselineClingoCmd <- paste(clingoCmd, clingoInputFiles, aspSetsFullPathDirect, ancestralResults)
      } else {
        #baselineClingoCmd <- paste(clingoCmd, clingoInputFiles, ancestralResults)
        baselineClingoCmd <- paste(clingoCmd, clingoInputFiles)
      }
        
      scoresG <- foreach (i = 1:length(G)) %dopar% {
        learn.asp.iterative.loop (n, graph_filename_template,  (i-1)%% n + 1 , (i-1)%/% n + 1, predicateToTest= "edge",
                                multipleMinimaTempDir, baselineClingoCmd, optimum_baseline=optimum_baseline, baselineTimeout=baselineTimeout, debugOutput=debugOutput, verbose=verbose)
      }
    }
  } else {
    scoresC <- list()
    for (i in 1:length(C)) {
      score <- learn.asp.iterative.loop (n, graph_filename_template,  (i-1)%% n + 1 , (i-1)%/% n + 1, predicateToTest= "causes",
                              multipleMinimaTempDir, baselineClingoCmd, optimum_baseline=optimum_baseline, debugOutput=debugOutput, verbose=verbose)
      if (returnNothingIfTimeout && score$timeout) {
        return (list(C=array(0, c(n,n)), G=array(0, c(n,n)), Ge=array(0, c(n,n)), Gs=array(0, c(n,n)), objective=0, timeout=TRUE))
      }
      scoresC[[length(scoresC)+1]] <- score
    }
    
    if (evalDirectCauses) {     
      ancestralResults <-  paste(currentDir, "/", graph_filename_template, "_ancestral.txt", sep="")
      sink(ancestralResults)
      for (i in 1:length(scoresC)){
        arg1 <- (i-1)%% n + 1 
        arg2 <- (i-1)%/% n + 1
        if(scoresC[[i]]$score > 0) {
          cat(":- not causes(",arg1,",", arg2,").\n", sep="")
        } else if(scoresC[[i]]$score < 0) {
          cat(":- causes(",arg1,",", arg2,").\n", sep="")
        }   
      }
      sink()
      
      if (grepl("tom", clingoInputFiles) > 0){
        clingoInputFiles <- gsub("ASP/tom_acyclic_complete.pl", "ASP/new_wmaxsat_acyclic_causes.pl", clingoInputFiles)
        clingoInputFiles <- gsub("ASP/tom_acyclic_maxschedule_complete.pl", "ASP/new_wmaxsat_acyclic_causes.pl", clingoInputFiles)
        aspSetsFullPathDirect <-  paste(currentDir, "/", graph_filename_template, "_direct.pre.asp", sep="")
        writeAspSets(n=n, aspSetsFullPath = aspSetsFullPathDirect, tested_independences=tested_independences)
        build_tree(n=n, aspSetsFullPath = aspSetsFullPathDirect, tested_independences=tested_independences)   
        baselineClingoCmd <- paste(clingoCmd, clingoInputFiles, aspSetsFullPathDirect, ancestralResults)
      } else {
        baselineClingoCmd <- paste(clingoCmd, clingoInputFiles, ancestralResults)
      }
      
      scoresG <- list()
      for (i in 1:length(G)) {
        score <- learn.asp.iterative.loop (n, graph_filename_template,  (i-1)%% n + 1 , (i-1)%/% n + 1, predicateToTest= "edge",
                                           multipleMinimaTempDir, baselineClingoCmd, optimum_baseline=optimum_baseline, debugOutput=debugOutput, verbose=verbose)
        if (returnNothingIfTimeout && score$timeout) {
          return (list(C=array(0, c(n,n)), G=array(0, c(n,n)), Ge=array(0, c(n,n)), Gs=array(0, c(n,n)), objective=0, timeout=TRUE))
        }
        scoresG[[length(scoresG)+1]] <- score
      }
    }
  }

  # The previous loop was in parallel, once it's done we can use the values to populate C.
  for (i in 1:length(C)){
    if (returnNothingIfTimeout && scoresC[[i]]$timeout) {
      # Unluckily cannot do in parallel foreach.
      return (list(C=array(0, c(n,n)), G=array(0, c(n,n)), Ge=array(0, c(n,n)), Gs=array(0, c(n,n)), objective=0, timeout=TRUE))
    }
    arg1 <- (i-1)%% n + 1 
    arg2 <- (i-1)%/% n + 1
    if (scoresC[[i]]$timeout) {
      C[arg2, arg1] <- 0
    } else {
      C[arg2, arg1] <- scoresC[[i]]$score
    }
    if (evalDirectCauses) {
      G[arg2, arg1] <- scoresG[[i]]$score
    } 
  }
  
  list(C=C, G=G, Ge=array(0, c(n,n)), Gs=array(0, c(n,n)), objective=optimum_baseline, timeout=FALSE)
}

learn.asp.iterative.loop <- function (n, graph_filename_template,
                                      arg1, arg2,
                                      multipleMinimaTempDir, originalClingoCmd,
                                      optimum_baseline,
                                      predicateToTest= "causes",
                                      baselineTimeout = FALSE,
                                      debugOutput=FALSE,
                                      returnNothingIfTimeout=FALSE,
                                      verbose=1) {
  if (arg1 == arg2){
    L <- list(score=-Inf, timeout=FALSE)
    return(L)
  }
  
  optimum_causes <- 0
  optimum_not_causes <- 0
  
  atLeastOneTimeout <- FALSE
  
  for (type in list("", "not")) {    
    optimum_or_timeout <- learn.asp.iterative.loop.internal(n, graph_filename_template, arg1, arg2,
                                      multipleMinimaTempDir, originalClingoCmd, predicateToTest=predicateToTest,
                                      type=type, debugOutput=debugOutput, verbose=verbose)
    
    if (optimum_or_timeout$timeout){
      atLeastOneTimeout <- TRUE
      if (returnNothingIfTimeout){
        return (list(score=0, timeout=TRUE))
      }
    }
    
    optimum <- optimum_or_timeout$optimum
    
    if (type != "not") {
      optimum_not_causes <- optimum
    } else {
      optimum_causes <- optimum
    }
    
    if (!baselineTimeout && optimum != optimum_baseline) {
      # No need to do the next minimization, we already found the score.
      next
    }
  }
  
  list(score=optimum_not_causes-optimum_causes, timeout=atLeastOneTimeout)
}

learn.asp.iterative.loop.internal <- function (n, graph_filename_template,
                                      arg1, arg2,
                                      multipleMinimaTempDir, originalClingoCmd,
                                      type,
                                      predicateToTest= "causes",
                                      debugOutput=FALSE,
                                      returnNothingIfTimeout=FALSE,
                                      verbose=1){
  
  extraAspBaseName <- paste(multipleMinimaTempDir, graph_filename_template, type, arg1, arg2, sep="_")
  extraAspProgram <- paste(extraAspBaseName, ".pl", sep="")
  cat(paste(":- ", type, " ", predicateToTest, "(", arg1, ",", arg2, ").",sep =""), file=extraAspProgram)
  
  if (debugOutput) {
    optCmd <- paste(originalClingoCmd, extraAspProgram, "|tee ", paste(extraAspBaseName, ".txt", sep=""))
  } else {
    optCmd <- paste(originalClingoCmd, extraAspProgram)
  }
  optResult <- system(optCmd, intern = TRUE)
  
  if ((length(grep("INTERRUPTED", optResult))> 0 || length(grep("TIME LIMIT", optResult)))> 0){
    # Timeout
    if (returnNothingIfTimeout) {return (list(optimum=0, timeout=TRUE))}
    optTimeout = TRUE
  } else {
    optTimeout = FALSE
  }

  if (length(grep("UNSATISFIABLE", optResult))> 0) {
    # The optimization is unsatisfiable.
    return (list(optimum=Inf, timeout=FALSE))
  } else if (length(grep("OPTIMUM FOUND", optResult))> 0){
    optimumLine <- optResult[grep('Optimization : ', optResult)]
    optimum <- as.numeric(gsub('Optimization : ','', optimumLine))
  } else if (length(grep("SATISFIABLE", optResult))> 0){
    optimumLine <- optResult[grep('Optimization : ', optResult)]
    optimum <- as.numeric(gsub('Optimization : ','', optimumLine))
  } else if (optTimeout) {
    optimum <- 0
  }  else {
    stop("Strange behaviour of solver:", optResult)
  }
  list(optimum=optimum, timeout=optTimeout)
}

#
#
#
########################## DEPRECATED ###################################

learn.asp.iterativeSAT <- function( n, J, graph_filename_template,
                                 multipleMinima = FALSE,
                                 solver_conf="--time-limit=25000 --quiet=1,0",
                                 currentDir = "./../tmp/",
                                 verbose=0,
                                 asp_program,
                                 aspSetsFullPath,
                                 outputFullPath,
                                 indFullPath,
                                 background_knowledge_file=NULL,
                                 debugOutput=FALSE){
  multipleMinimaTempDir <- paste(currentDir, "/multipleMinimaTemp/", sep ="")
  system(paste("mkdir -p", multipleMinimaTempDir))
  multipleMinimaIterativeC <- array(2, c(n,n))
  
  satAspProgram <- paste(multipleMinimaTempDir, graph_filename_template, sep="")
  system(paste("sed -e 's/.*minimize.*//g' -e 's/fail([^)]*)//g'", asp_program, "> ", satAspProgram))
  
  for (arg1 in 1:n) {
    for (arg2 in 1:n) {
      if (arg1 == arg2){
        multipleMinimaIterativeC[arg1, arg2] <- 0
        next
      }
      for (type in list(":-", "")) {
        extraAspProgram <- paste(paste(multipleMinimaTempDir, graph_filename_template, type, arg1, arg2, sep="_"), ".pl", sep="")
        # The following is used for the oracle version: 
        #system(paste("echo '\n", type, "causes(", arg1, ",", arg2, ").", "\n' >> ", extraAspProgram, sep =""))
        cat(paste(type, "causes(", arg1, ",", arg2, ").",sep =""), file=extraAspProgram)
        # "clingo" was "./../ASP/clingo430"
        clingoCmd <- paste("clingo", solver_conf, aspSetsFullPath, indFullPath, satAspProgram, extraAspProgram)
        if (debugOutput) {
          extraOutputFullPath <- paste(paste(multipleMinimaTempDir, graph_filename_template, type, arg1, arg2, sep="_"), ".ind.clingo", sep="")
          unSATcmd <- paste(clingoCmd, "| tee" , extraOutputFullPath, "|grep UNSATISFIABLE")
        } else {
          unSATcmd <- paste(clingoCmd, "|grep UNSATISFIABLE")
        }
        if (verbose > 10 ) cat(clingoCmd, "\n")
        # Grep returns 1 if there is no UNSATISFIABLE, i.e. if it is satisfiable:
        satisfiable <- system(unSATcmd)
        if (!satisfiable) {
          # If type is "" there is no model when causes(arg1,arg2) is true => causes is always false.
          # If type is ":-" there is no model when :-causes(arsg1,arg2) is true => causes is (always or sometimes??) true.
          multipleMinimaIterativeC[arg2, arg1] <- if (type ==":-") {1} else {0}
        }
      }
    }     
  }
  if (verbose) {
    cat("\nLearned C\n")
    print(multipleMinimaIterativeC)
  }
  L <- list()
  L$C <- multipleMinimaIterativeC
  L$G <- multipleMinimaIterativeC
  L$Ge <- array(0, c(n,n))
  L$Gs <- array(0, c(n,n))
  L$objective <- 0
  L$compareC <- C
  L
}

