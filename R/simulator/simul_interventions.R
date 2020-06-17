simul_interventions.generateModel <- function(p=6, numConf=0, numInts = 0,
                                             # Each value is the probability of the number of children, e.g. 0.5 is the prob of 2 children
                                             # In this case there cannot be more than 4 children.
                                             prob=c(0.125,0.5,0.25,0.125), 
                                             probConf=c(0,0,1),
                                             probInts=c(0.125,0.5,0.25,0.125),
                                             weights=0,
                                             debug=T) {
  # Total number of variables (regime + confounders + intervention variables + observed )
  P <- 1 + numConf + numInts + p
  
  # First, a random directed acyclic graph representing the true causal structure is built.
  # B[i,j] = coefficient of linear combination.
  B <- matrix(0,P,P)
  
  # We start with numConf "global" confounders, i.e., unobserved variables that influence observed ones.
  if (numConf > 0 ) {
    for (i in 1:numConf) {
      Nchildren <- sample(x=length(probConf),size=1, prob=probConf)
      Nchildren <- min(p, Nchildren-1)
      cat('Confounder', i+1, 'has', Nchildren, 'children:\n')
      if( Nchildren > 0 ) {
        ch <- sample(x=p, size=Nchildren,  replace=FALSE)
        if (debug) cat('- children of', i+1, '= ', paste(numConf + numInts + ch + 1, collapse=","), '\n')
        randomweights <- list()
        if( weights == 1 ) {
            randomweights <- runif(Nchildren, -sqrt(3), sqrt(3))
          } else if( weights == -1 ) {
            randomweights <- runif(Nchildren, 2./3., 4./3.) * (rbinom(Nchildren, 1, 0.5)*2-1)
          } else if( weights == -2 ) {
            randomweights <- runif(Nchildren, 0, 2*sqrt(3))
          } else {
            randomweights <- rnorm(Nchildren, 0.2,0.8) * sample(c(-1,1),1)
          }
        B[i + 1, numConf + numInts + ch + 1] <- randomweights
      }
    }
  }
  
  # Similarly intervention variables cause observed ones, but never the opposite.
  if (numInts > 0 ) {
    for (i in 1:numInts) {
      Nchildren <- sample(x=length(probInts), size=1, prob=probInts)
      Nchildren <- min(p, Nchildren-1)
      cat('Intervention var', numConf + i + 1, 'has', Nchildren, 'children\n')
      if ( Nchildren > 0 ) {
        ch <- sample(x=p, size=Nchildren,  replace=FALSE)
        if (debug) cat('- children of', numConf + i + 1, '= ', paste(numConf + numInts + ch + 1, collapse=","), '\n')
        randomweights <- list()
          if( weights == 1 ) {
            randomweights <- runif(Nchildren, -sqrt(3), sqrt(3))
          } else if( weights == -1 ) {
            randomweights <- runif(Nchildren, 2./3., 4./3.) * (rbinom(Nchildren, 1, 0.5)*2-1)
          } else if( weights == -2 ) {
            randomweights <- runif(Nchildren, 0, 2*sqrt(3))
          } else {
            randomweights <- rnorm(Nchildren, 0.2,0.8) * sample(c(-1,1),1)
          }       
        B[numConf + i + 1, numConf + numInts + ch + 1] <- randomweights
      }
      # Each intervention is a deterministic function of the regime.
      B[1, numConf + i + 1] <- rnorm(1)
    }
  }
  
  for (i in p:1) {
    if (all(prob == 0)) {
      Nchildren <- 0
    } else {
      Nchildren <- sample(x=length(prob), size=1, prob=prob)
      Nchildren <- min(p, Nchildren-1)
    }
    if (Nchildren > (p - i)) {
        Nchildren <- p - i
    }
    if (debug) cat('Standard var', numConf + numInts + i + 1 , 'has', Nchildren, 'children\n')
    if (Nchildren > 0) {
      ch <- p + 1 - sample(x = p - i, size=Nchildren, replace=FALSE)
      if (debug) cat('- children of', numConf + numInts + i + 1, '= ', paste(numConf + numInts + ch + 1, collapse=","), '\n')
      randomweights <- list()
      if( weights == 1 ) {
          randomweights <- runif(Nchildren, -sqrt(3), sqrt(3))
        } else if( weights == -1 ) {
          randomweights <- runif(Nchildren, 2./3., 4./3.) * (rbinom(Nchildren, 1, 0.5)*2-1)
        } else if( weights == -2 ) {
          randomweights <- runif(Nchildren, 0, 2*sqrt(3))
        } else {
          randomweights <- rnorm(Nchildren, 0.2,0.8) * sample(c(-1,1),1)
      }
      B[numConf + numInts + i + 1, numConf + numInts + ch + 1] <- randomweights
    }
  }
  
  # Apply random permutation of variables, if requested
  perm <- sample(1:p)

  permutedB <- B[(numConf+numInts+2):P,(numConf+numInts+2):P]
  permutedB <- permutedB[perm,perm]

  ### Convert to standard format data.
  
  gen_M <- list(G=array(0,c(p,p)), Ge=array(0,c(p,p)), Gs=array(0,c(p,p)))
  gen_M$perm <- perm
  gen_M$B <- t(permutedB)
  
  fullperm <- c((numConf+numInts+1+perm), 1, (numConf+2):(numConf+numInts+1))
  gen_M$fullB <- t(B[fullperm,fullperm])
  gen_M$fullGe <- array(0,c(nrow(gen_M$fullB),nrow(gen_M$fullB)))
  gen_M$fullperm <- fullperm
  
  gen_M$G <- gen_M$B
  gen_M$G[which(gen_M$G != 0)] <- 1
  gen_M$C <- dag2causes(gen_M$G)
  
  if (numConf > 0) {
    for (i in 1:numConf) {
      childrenOfI <- B[i + 1, (numConf+numInts+2):P]
      confoundedB <- childrenOfI[perm]
      confPart <- which(confoundedB != 0)
      if (length(confPart)> 1){
        confPairs <- t(combn(confPart, 2))
        for (j in 1:nrow(confPairs)) {
          gen_M$Ge[confPairs[j,1], confPairs[j,2]] <- 1
          gen_M$Ge[confPairs[j,2], confPairs[j,1]] <- 1
        }  
      }
      childrenOfI <- B[i + 1, ]
      fullconfoundedB <- childrenOfI[fullperm]
      confPart <- which(fullconfoundedB != 0)
      if (length(confPart)> 1){
        confPairs <- t(combn(confPart, 2))
        for (j in 1:nrow(confPairs)) {
          gen_M$fullGe[confPairs[j,1], confPairs[j,2]] <- 1
          gen_M$fullGe[confPairs[j,2], confPairs[j,1]] <- 1
        }  
      }
    }
  }
  
  gen_M$fullG <- gen_M$fullB
  gen_M$fullG[which(gen_M$fullG != 0)] <- 1
  gen_M$fullC <- dag2causes(gen_M$fullG)
  gen_M$fullG[which(gen_M$fullG == 0)] <- -1
  gen_M$fullC[which(gen_M$fullC == 0)] <- -1

  gen_M$original_B <- B
  gen_M$numConf <- numConf
  gen_M$numInts <- numInts

  gen_M
}

simul_interventions.generateData <- function(M, disturbance_distribution=function(x){0.1*rnorm(x,0.2,0.8)}, base_function=identity, 
                                            nObs=500, weights=0) {
  B <- M$original_B
  numConf <- M$numConf
  numInts <- M$numInts
  p <- nrow(M$C)
  
  P <- 1 + numConf + numInts + p

  ## Simulate observational data
  ##allObs <- matrix(rnorm(nObs*P),nrow=nObs) %*% diag(scales * noiseSigma)
  ##X <- allObs / (eye(P) - B)   # SLOW
  #obs <- t(forwardsolve(t(diag(1,P) - B), t(allObs), transpose=FALSE))
  
  # Simulate observational data
  alldata <- matrix(0, P, nrow=nObs)
  # Not very efficient: one by one.
  for (j in 1:nObs) {
    V <- array(0, c(P))
    V[1] <- sample(x = numInts + 1, 1)
    for (i in 2:P) {
      if ( numInts > 0 && i > numConf + 1 && i <= numConf + numInts + 1 ) {
        # Intervention variables are deterministic functions of the regime.
        # For now diagonal design:
        if (V[1] == i - numConf) {
          value_i <- 1
        } else {
          value_i <- 0
        }
      } else {
        disturbance_i <- disturbance_distribution(1);
        value_i <- disturbance_i;
        # If it has parents, count their contribution.
        parents <- which(B[,i] != 0)
        for (parent in parents) {
            value_i <- value_i + B[parent,i] * base_function(V[parent])
        }
      }
      V[i] <- value_i
    }
    alldata[j, ] <- V
  }
 
  # Only select observed variables
  obs <- alldata[,(numConf+numInts+2):P]
  
  # Apply random permutation of variables, if requested
  perm <- M$perm
  obs <- obs[,perm]

  regime <- alldata[,1]
  ints <- alldata[,(numConf+2):(numConf+numInts+1)]
  
  final_data <- data.frame(cbind(obs, regime))
  names(final_data) <- c(paste("v", 1:p, sep=""), "experiment")
  
  ### Convert to standard format data
  D <- list()
  E <- experimentConfiguration(p, "passive")
  D[[1]] <- list()
  D[[1]]$M<-M
  D[[1]]$e<-E[1,]
  D[[1]]$N<-nObs
  D[[1]]$all_data <- alldata
  D[[1]]$data <- final_data
  D[[1]]$Cx<-cov(D[[1]]$data)
  
  
  list(M=M, D=D)
}
