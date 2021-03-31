pipeline<-function(n=4, topology="random", # model properties
                   exconf="passive", N=500, # sampling properties
                   pedge = 1/((n-1)*2), # probablity of edge in the generated model.
                   restrict = c('acyclic'),
                   confounder_proportion=0.5,
                   test="bayes", schedule=1, # testing properties
                   p=0.5, alpha=1.5, # eq. sample size for the bayes test
                   weight="log", 
                   encode="new_wmaxsat_acyclic_causes.pl", # how to encode them
                   #encode=NULL,
                   solver="clingo", # which solver to use
                   solver_conf="--time-limit=1000 --quiet=1",
                   #solver_conf="--configuration=crafty --time-limit=25000 --quiet=1,0",
                   multipleMinima = "iterative",
                   model=NULL, # input model if given
                   samples=NULL, # input sample data if given
                   indPath = NULL, # input .ind file if given
                   background_knowledge_file=NULL, # background knowledge file if available.
                   # what restrictions apply to model space, note that encode file 
                   # tells which restrictions are forced by the learned model space,
                   intervened_variables = c(),
                   tmpDir = "../tmp/",
                   outputCsv=TRUE,
                   outputClingo = NULL,                   verbose=1){ # how much information to print
  #Pipeline for running the algorithm inference. This first creates model and data, then runs
  #the requested inference and finally assesses and prints out the quality of the inference.
  #Use "set.seed()" before this function to compare the performance of several algorithms.
  #
  # Model properties:
  ###################
  # n  - number of observed variables (default=4)
  #
  #Data properties:
  #################
  #exconf   - experiment configuration, "passive" (default) or "single" or "random"
  # N       -total number of data samples (divided equally to all experiments) (default=1000)
  #
  #Independence test properties:
  ##############################
  #test     -type of test used to generate data, outputs a value between 0 and 1?
  #         -"classic" classic correlation test (p-value = 0.05)
  #         -"oracle" independence facts determined by Patrik's oracle
  #         -"bayes" integrating over the parameters
  #         -"BIC" (default) BIC-based approximation of the bayes test, faster and almost as good
  #         -the prior parameters are a little different, only p as the prior prob. is needed.
  #schedule -independence test scheduling
  #         -maximum conditioning set size for the independence tests
  #         -if a number T it means all tests up to intervention tests size T
  #           (can be Inf, default n-2 so all possible conditioning test sizes)
  #weight   -how to determine the weights from the dependencies:
  #         -"log" take log of the probability
  #         -"constant" the most likely option (dep or indep) gets weight 1
  #         -"hard-deps" put dependencies as hard constraints
  #         -for competing algorithms, pc and so on, use any of the above.
  ###############################
  #encode   - which encoding to use, gives the file in the ASP/ directory
  #solver   -"clingo", or "pcalg-fci","pcalg-cfci" or "bnlearn" for score based learning
  #p        -for bayes and BIC tests the apriori probability of 
  #         -for classic test using algorithms the p-value threshold
  #         -for BIC-based score based learning the prior parameter
  #alpha    -for bayes test the eq. sample size
  ###############################
  # solver_conf  - a string which defines additional parameters to clingo
  #Printing options:
  #verbose  -0 to not print anything, 1 to some printing
  ##############################################################################
  
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
  testName <- paste(test, paste(intervened_variables, collapse="",sep=""), sep='')
  graph_filename_template <-paste(n, topology, N, testName, "sched", paste(schedule, collapse=""), weight, gsub("./","", encode) , solver, "multipleMin", multipleMinima, p, alpha, sep="_")
  
  # Simulate or use the provided data.
  MD <- simulateData(simulationConfig=simulationConfig, samples=samples, model=model, indPath=indPath, 
                     isOracle= (test=="oracle"), verbose=verbose) 
  
  # Perform independence tests or read them from .ind.
  # Note: if we use a solver from the pcalg package the tests are supposed to be done incrementally by the algorithm, so we will fake it by using
  # the tests we performed here.
  if (!is.null(indPath)) {
    tic()
    parsed_indeps <- parse_asp_indeps(indPath)
    tested_independences <- parsed_indeps$tested_independences
    testConfig$n <- parsed_indeps$n
    if (is.null(testConfig$n)) {testConfig$n <- n}
    test_time <- toc()
  } else {
    tic()
    tested_independences_and_indpath <- conductIndepTests(D=MD$D, testConfig=testConfig, verbose=verbose)
    tested_independences <- tested_independences_and_indpath$tested_independences
    indPath <- tested_independences_and_indpath$indPath
    test_time <- toc()
    
  }
  
  if (is.null(encode)) { 
    # Use only to compute tested independences.
    return (list(solving_time=0, testing_time=test_time, encoding_time=0, objective=NA, tested_independences=tested_independences))
  }
  
  # Try to learn a model.
  evalDirectCauses = TRUE
  L <- learn(testConfig, solverConfig, currentDir = currentFolder, graph_filename_template=graph_filename_template, indPath=indPath, 
             background_knowledge_file=background_knowledge_file, tested_independences=tested_independences, MD=MD,
             evalDirectCauses=evalDirectCauses, verbose=verbose)
  L$tested_independences <- tested_independences
  L$testing_time <- test_time
  # Save also the true model.
  L$M <- MD$M
  
  # If the solving time is infinite, the model is not good.
  if (is.infinite(L$solving_time)) {
    return(L)
  }
  
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
    write.csv(x=MD$D[[1]]$data, file=paste(currentFolder, graph_filename_template, "_data.csv", sep=""))
  }
  
  # Works also for FCI.
  if (!is.null(outputClingo)){
    G <- if (evalDirectCauses) {L$G} else {L$C} 

    sink(outputClingo)
    for (i in 1:nrow(G)){
      for (j in 1:ncol(G)){
        if (G[i,j] < 0) {
          cat("-")
        } 
        cat("causes(",j,",", i,")=", abs(G[i,j]), sep="")
        if (!is.infinite(G[i,j])) cat(".00000", sep="")
        cat("\n", sep="")
      }
    }
    sink()
  }
  
  invisible(L)
}

