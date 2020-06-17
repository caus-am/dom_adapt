simulateData <- function(simulationConfig, samples, model, indPath, isOracle = FALSE, typeOfSimulator="random", returnData=TRUE, verbose=0) {
  # Function that simulates the model and then samples the data, if necessary.
  # Returns a list with M=true model, D=sampled data.
  # - simulationConfig$n  
  # - simulationConfig$topology
  # - simulationConfig$exconf
  # - simulationConfig$N
  # - simulationConfig$pedge
  # - simulationConfig$restrict
  # - samples: NULL unless the data is already provided.
  # - model: NULL unless the true model is already provided.
  # - indPath: NULL unless the independence tests have already been ran. In that case, we don't need to generate anything.
  # - isOracle: FALSE for noisy data, TRUE for oracle data. 
  
  # Independence test are already done, there is no need to simulate anything.
  if (!is.null(indPath)) {
    return (list(M=model, D=NULL))
  }
  
  # No model or sampling necessary.
  # Assume only observational data (needs revision otherwise).
  if (!is.null(samples)) {
    D <- list()
    D[[1]] <- list()
    D[[1]]$M <- model
    D[[1]]$data <- samples
    D[[1]]$e <- rep(0, ncol(samples))
    D[[1]]$N <- nrow(samples)
    return (list(M=model, D=D))
  }
  
  # Get the true model, either from model or by generating a new one randomly.
  if (!is.null(model) ) {
    M <- model
  } else {
    cat("\n* Generating the model: n=", simulationConfig$n,", restrict=", simulationConfig$restrict, 
          ", topology=",simulationConfig$topology, ", probability of edge = ", simulationConfig$pedge, 
          ", confounder proportion=",simulationConfig$confounder_proportion,
          ".\n",sep='')
    
    if (typeOfSimulator == "int") {
      M <- simul_interventions.generateModel(p=simulationConfig$n,
                                             numConf = as.integer(simulationConfig$n*simulationConfig$confounder_proportion),
                                             numInts = simulationConfig$numInts)
    } else if (typeOfSimulator == "int_example") {
      M <- simul_interventions.generateModel(p=2, prob=c(0, 0, 0, 0), probInts = c(0, 1),
                                             numConf = 0,
                                             numInts = 1)
    } else if (typeOfSimulator == "int_example2") {
      M <- simul_interventions.generateModel(p=5, prob=c(0.25, 0.25, 0.25, 0.25), probInts = c(0.25, 0.5, 0.25),
                                             numConf = 1,
                                             numInts = 2)
    } else {
      M <- generateModel(n=simulationConfig$n, restrict=simulationConfig$restrict, topology=simulationConfig$topology, 
                       model=model, samples=samples, pedge=simulationConfig$pedge, confounder_proportion=simulationConfig$confounder_proportion, 
                       verbose=verbose)
    }
  }
  
  if (returnData) {
    if (isOracle) {
      if (verbose) cat(" - Skipping sample data generation, since using oracle.\n",sep='')  
      D <- list()
      E <- experimentConfiguration(simulationConfig$n, simulationConfig$exconf)
      # For each experiment vector e in E:
      for ( i in 1:nrow(E)) {
        # Create the tuple that will be stored in D, storing the vector e.
        D[[i]]<-list(e=E[i,],M=M)
        # Get indexes of intervened variables in this particular experimental setting.
        J<-which( D[[i]]$e==1 )
        # The data consist of manipulated graphs where
        # Edge heads into the intervened variables are cut.
        D[[i]]$M$G[J,]<-0
        D[[i]]$M$Ge[J,]<-0
        D[[i]]$M$Ge[,J]<-0
      }
    } else if (grepl("int", typeOfSimulator)) {
      MD <- simul_interventions.generateData(M=M, nObs=simulationConfig$N)
      D <- MD$D
    } else {
      if (verbose) {
        cat("\n* Generating the sample data using the experiment configuration: exconf=", simulationConfig$exconf, 
            ", number of samples=", simulationConfig$N, ".\n",sep='')
      }
      D <- generateSampleData(n=simulationConfig$n, exconf=simulationConfig$exconf, test=test, N=simulationConfig$N, 
                              M=M, verbose=verbose)
    }
  } else {
    D <- list()
  }
  # Change for compatibility with results from the learn step: not causes = -1.
  M$G[(M$G==0)] <- -1
  M$Ge[(M$Ge==0)] <- -1
  M$Gs[(M$Gs==0)] <- -1
  M$C[(M$C==0)] <- -1
  
  list(M=M, D=D)
}

simulateData._test1 <- function() {
  simulationConfig <- list(n=4, topology="random", exconf="passive", N=500, pedge = 1/3, restrict = c('acyclic'))
  simulateData(simulationConfig=simulationConfig, samples=NULL, model=NULL, indPath=NULL, 
                     isOracle=0, verbose=0) 
}

simulateData._test2 <- function() {
  simulationConfig <- list(n=4, topology="random", exconf="single", N=500, pedge = 1/3, restrict = c('acyclic'))
  simulateData(simulationConfig=simulationConfig, samples=NULL, model=NULL, indPath=NULL, 
               isOracle=0, verbose=0) 
}

simulateData._test3 <- function() {
  simulationConfig <- list(n=4, topology="random", exconf="single", N=500, pedge = 1/3, restrict = c('acyclic'))
  simulateData(simulationConfig=simulationConfig, samples=NULL, model=NULL, indPath=NULL, 
               isOracle=1, verbose=0) 
}