simulateData_NIPS <- function(simulationConfig, samples, model, indPath, isOracle = FALSE, typeOfSimulator="random", filter_require_path_IY=FALSE, IFactor=1.0, returnData=TRUE, verbose=0) {
  # Function that simulates the model and then samples the data, if necessary.
  # Returns a list with M=true model, D=sampled data.
  # - simulationConfig$n  
  # - simulationConfig$topology
  # - simulationConfig$exconf
  # - simulationConfig$N
  #   OR: $NObs, NPerInt
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
    #browser()
    if (verbose) {
      cat("\n* Generating the model: n=", simulationConfig$n,", restrict=", simulationConfig$restrict, 
            ", topology=",simulationConfig$topology, ", probability of edge = ", simulationConfig$pedge, 
            ", confounder proportion=",simulationConfig$confounder_proportion,
            ".\n",sep='')
    }
    
    if (typeOfSimulator == "hughes") {
      M <- simulHughesData.generateModel('',topology=simulationConfig$topology, p=simulationConfig$n,
                                         nObs=simulationConfig$N,nInt=0,numConf=as.integer(simulationConfig$n*simulationConfig$confounder_proportion))
    } else if (typeOfSimulator == "exp") {
      M <- simul_experimental.generateModel(topology=simulationConfig$topology, p=simulationConfig$n,nObs=simulationConfig$N,numConf=as.integer(simulationConfig$n*simulationConfig$confounder_proportion))
    } else if (typeOfSimulator == "int") {
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
    } else if (typeOfSimulator == "int_NIPS2017") {
      probNumConf <- simulationConfig$probNumConf
      numConf <- sample(x=length(probNumConf), size=1, prob=probNumConf) - 1
      # The following parameters could be added in the call below to tweak graph generation (numbers of children per node):
      # - prob
      # - probConf
      # - probInts
      M <- simul_interventions_NIPS.generateModel(p=simulationConfig$n,
                                             numConf = numConf,
                                             numInts = simulationConfig$numInts,
                                             blind=TRUE,
                                             filter_require_path_IY=filter_require_path_IY,
                                             IFactor=IFactor,
                                             debug=verbose)
    } else if (typeOfSimulator == "int_NIPS2017_FSfails_5node") {
      # Example for which feature selection fails, and JCI can in principle detect the problem (from example2.tex)
      B <- matrix(0,6,6)
      # Arrows from R to I's
      B[1,2] <- B[1,3] <- 1
      # Arrows from I's to X's
      B[2,6] <- B[3,4] <- 1
      # Arrows between X's
      B[4,5] <- B[5,6] <- 1
      M <- simul_interventions_NIPS.generateModel(B=B, p=3, numConf=0, numInts=2,
                                                  blind=TRUE, IBlind=1, YBlind=2, shuffleObs=FALSE,
                                                  debug=verbose)
    } else if (typeOfSimulator == "int_NIPS2017_doubleLCD" || typeOfSimulator == "int_NIPS2017_FSfails_4node") { # 2nd name for compatibility
      # Double LCD example (from example.tex), where JCI should detect something (but we don't know what FS will do)
      B <- matrix(0,5,5)
      # Arrows from R to I's
      B[1,2] <- B[1,3] <- 1
      # Arrows from I's to X's
      B[2,4] <- B[3,4] <- 1
      # Arrows between X's
      B[4,5] <- 1
      M <- simul_interventions_NIPS.generateModel(B=B, p=2, numConf=0, numInts=2,
                                                  blind=TRUE, IBlind=1, YBlind=2, shuffleObs=FALSE,
                                                  debug=verbose)
    } else if (typeOfSimulator == "int_NIPS2017_FSfails_4node_Thijs") {
      # Example for which feature selection fails, and JCI can't detect the problem (Thijs' example)
      B <- matrix(0,5,5)
      # Arrows from R to I's
      B[1,2] <- 1
      # Arrows from I's to X's
      B[2,3] <- B[2,5] <- 1
      # Arrows between X's
      B[3,4] <- B[4,5] <- 1
      M <- simul_interventions_NIPS.generateModel(B=B, p=3, numConf=0, numInts=1,
                                                  blind=TRUE, IBlind=1, YBlind=2, shuffleObs=FALSE,
                                                  debug=verbose)
    } else if (typeOfSimulator == "int_NIPS2017_FSfails_4node_Thijs_scatterplot") {
      # Same example as above, with fixed values of B to get the scatterplots for the introduction
      B <- matrix(0,5,5)
      # Arrows from R to I's
      B[1,2] <- 1
      # Arrows from I's to X's
      B[2,3] <- .2 #.1
      B[2,5] <- .45 #.6
      # Arrows between X's
      B[3,4] <- .8 #1 #.8
      B[4,5] <- 1.6 #1.4
      M <- simul_interventions_NIPS.generateModel(B=B, p=3, numConf=0, numInts=1,
                                                  weights=NULL,
                                                  blind=TRUE, IBlind=1, YBlind=2, shuffleObs=FALSE,
                                                  debug=verbose)
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
    } else if (typeOfSimulator == "hughes") {
      MD <- simulHughesData.generateData('',M=M, p=simulationConfig$n,nObs=simulationConfig$N,nInt=0)
      D <- MD$D
    } else if (typeOfSimulator == "exp") {
      MD <- simul_experimental.generateData(M=M, p=simulationConfig$n,nObs=simulationConfig$N)
      D <- MD$D
    } else if (grepl("int", typeOfSimulator)) {
      MD <- simul_interventions_NIPS.generateData(M=M, nObs=simulationConfig$NObs, nPerInt=simulationConfig$NPerInt)
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
  #M$C[(M$C==0)] <- -1
  
  list(M=M, D=D)
}


simulateData.NIPS2017 <- function(n=3, numInts=2, NObs=1000, NPerInt=1000, typeOfSimulator="int_NIPS2017", filter_require_path_IY=FALSE, IFactor=1.0) {
  # Generate mouse-challenge-type data
  # * Needs source("simul_interventions_NIPS.R")
  # * New values in output:
  #   M$IBlind: intervention being blinded, counting from 1 for first interventional regime (R=2);
  #   M$YBlind: observed variable being blinded, counting from 1 for first observed var in permuted data.
  #   D[[1]]$data_blinded: like D[[1]]$data, but with nan's where values need to be predicted.

  # Other choices for typeOfSimulator ("int_NIPS2017_*") give fixed example graphs (but still random weights and data)

  simulationConfig <- list(n=n, numInts=numInts, probNumConf=c(.5,.25,.25), NObs=NObs, NPerInt=NPerInt)
  simulateData_NIPS(simulationConfig=simulationConfig, samples=NULL, model=NULL, indPath=NULL, 
               isOracle=0, typeOfSimulator=typeOfSimulator, filter_require_path_IY=filter_require_path_IY, IFactor=IFactor,
               verbose=1)
}

simulateData.NIPS2017.large <- function(typeOfSimulator="int_NIPS2017", filter_require_path_IY=FALSE, IFactor=1.0) {
  simulateData.NIPS2017(n=4, numInts=2, NObs=NObs, NPerInt=NPerInt,
                        typeOfSimulator=typeOfSimulator, filter_require_path_IY=filter_require_path_IY, IFactor=IFactor)
}
