bootstrap<-function(repeat_bootstrap=10, n=4, topology="random", # model properties
                   exconf="passive", N=500, # sampling properties
                   pedge = 1/(n-1), # probablity of edge in the generated model.
                   confounder_proportion=0.5,
                   restrict = c('acyclic'),
                   test="bayes", schedule=n-2, # testing properties
                   p=0.5, alpha=1.5, # eq. sample size for the bayes test
                   weight="log", 
                   encode="new_wmaxsat_acyclic_causes.pl", # how to encode them
                   solver="clingo", # which solver to use
                   solver_conf="--time-limit=1000 --quiet=1",
                   intervened_variables = c(),
                   #solver_conf="--configuration=crafty --time-limit=25000 --quiet=1,0",
                   multipleMinima = "iterative",
                   model=NULL, # input model if given
                   samples=NULL, # input sample data if given
                   tmpDir = "../tmp/",
                   outputCsv=TRUE,
                   verbose=1){ # how much information to print
  # Same interface as pipeline.R
  
  if (tmpDir == "../tmp/") {
    currentFolder <- paste(tmpDir, as.numeric(as.POSIXct( Sys.time())), "/", sep="")
  } else {
    currentFolder <- paste(tmpDir, "/", sep="")
  }
  simulationConfig <- list(n=n, topology=topology, exconf=exconf, N=N, pedge = pedge, restrict = restrict, confounder_proportion=confounder_proportion)
  
  if (is.na(p)) p <- 0.5
  testConfig <- list(n=n, test=test, schedule=schedule, p=p, alpha=alpha, weight=weight, conditioning_vars=NULL, currentDir=currentFolder)
  solverConfig <- list(solver=solver, encode=encode, solver_conf=solver_conf, multipleMinima=multipleMinima, 
                       intervened_variables=intervened_variables)
  
  system(paste("mkdir -p", currentFolder))
  # Use a template for all the files related to each run.
  testName <- paste(test, paste(c(), collapse="",sep=""), sep='')
  graph_filename_template <-paste(n, topology, N, testName, "sched", paste(schedule, collapse=""), weight, gsub("./","", encode) , solver, "multipleMin", multipleMinima, p, alpha, sep="_")
  
  # Simulate or use the provided data.
  MD <- simulateData(simulationConfig=simulationConfig, samples=samples, model=model, indPath=NULL, 
                     isOracle= (test=="oracle"), verbose=verbose) 
  
  L <- list()
  # Save also the true model.
  L$M <- MD$M
  L$C <- array(0, c(n,n))
  L$G <- array(0, c(n,n))
  L$Ge <- array(0, c(n,n))
  L$solving_time <- 0
  
  # Perform independence tests.
  #for(bootstrap_iter in 1:repeat_bootstrap) {
  foreach (bootstrap_iter=1:repeat_bootstrap) %do% {
    # Pick randomly half of the samples.
    indices <- as.integer(runif(N/2)*N/2)
    halfD <- MD$D
    halfD[[1]]$data <- MD$D[[1]]$data[indices, ]
    halfD[[1]]$N <- halfD[[1]]$N/2
    #halfD[[1]]$Cx <- cov(halfD[[1]]$data)
    halfMD <- list(D=halfD)
    
    tested_independences_and_indpath <- conductIndepTests(D=halfD, testConfig=testConfig, verbose=verbose)
    tested_independences <- tested_independences_and_indpath$tested_independences
    indPath <- tested_independences_and_indpath$indPath
  
    # Try to learn a model.
    P <- learn(testConfig, solverConfig, currentDir = currentFolder, graph_filename_template=graph_filename_template, indPath=indPath, 
             background_knowledge_file=NULL, tested_independences=tested_independences, MD=halfMD, verbose=verbose)
    L$C <- L$C + P$C
    L$G <- L$G + P$G
    L$Ge <- L$Ge + P$Ge
    L$solving_time <- L$solving_time + P$solving_time
  }
  L$objective <- -1;
  
  L$C <- L$C/repeat_bootstrap
  L$G <- L$G/repeat_bootstrap
  L$Ge <- L$Ge/repeat_bootstrap
  
  # If the solving time is infinite, the model is not good.
  if (is.infinite(L$solving_time)) {
    return(L)
  }
  
  # Print the learned and true graph to dot files.
  plot_graph_as_dot(learnedGraph=L, trueGraph=MD$M, graph_filename_template=graph_filename_template, 
                    dotFilesFolder=currentFolder, names=colnames(MD$D[[1]]$data), verbose=verbose)
  
  if (verbose) {
    cat("Printing true and learned model.")
    if (!is.null(MD$M)){
      cat("\nTrue model G:\n")
      print(MD$M$G)
      cat("\nTrue model C:\n")
      print(MD$M$C)
    }
    cat("Learned G:\n")
    print(L$G)
    cat("Learned Ge:\n")
    print(L$Ge)
    cat("Learned C:\n")
    print(L$C)
  }
  if (outputCsv) {
    write.csv(x=L$C, file=paste(currentFolder, graph_filename_template, ".csv", sep=""))
    write.csv(x=L$G, file=paste(currentFolder, graph_filename_template, "_direct.csv", sep=""))
  }
  
  invisible(L)
}

