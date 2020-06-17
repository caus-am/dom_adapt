simul_interventions_NIPS.generateB <- function(P, p, numConf, numInts, prob, probConf, probInts, weights, debug) {
  # First, a random directed acyclic graph representing the true causal structure is built.
  # B[i,j] = coefficient of linear combination.
  B <- matrix(0,P,P)
  
  # We start with numConf "global" confounders, i.e., unobserved variables that influence observed ones.
  if (numConf > 0 ) {
    for (i in 1:numConf) {
      Nchildren <- sample(x=length(probConf),size=1, prob=probConf)
      Nchildren <- min(p, Nchildren-1)
      if (debug) cat('Confounder', i+1, 'has', Nchildren, 'children:\n')
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
      # I never allow an intervention variable to affect *all* observed variables
      Nchildren <- min(p-1, Nchildren-1)
      if (debug) cat('Intervention var', numConf + i + 1, 'has', Nchildren, 'children\n')
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
      B[1, numConf + i + 1] <- 1 # placeholder value: not a deterministic function
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
    if (debug > 2) cat('Standard var', numConf + numInts + i + 1 , 'has', Nchildren, 'children\n')
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

  B
}

# TODO: transpose is probably wrong
# dag2causes <- function( G ) {
#   # Calculates the causes-matrix for a DAG
#   # tc[i,j] == 0 if i is not a causal ancestor of j
#   # tc[i,j] == 1 if i is a causal ancestor of j
#   # (this is basically the transitive closure, without the diagonal)
# 
#   # Originally from graph_utils; rewritten to not rely on RBGL
#   n <- nrow(G)
#   tc <- G
#   # Floyd-Warshal algorithm for transitive closure
#   for (k in 1:n) {
#     for (i in 1:n) {
#       for (j in 1:n) {
#         if (tc[i,k] && tc[k,j]) {
#           tc[i,j] <- 1
#         }
#       }
#     }
#   }
#   tc <- t(tc)
#   diag(tc) <- 0
#   tc
# }

simul_interventions_NIPS.generateModel_unfiltered <- function(B=NULL, p=4, numConf=2, numInts = 2,
                                             # Each value is the probability of the number of children, e.g. 0.5 is the prob of 1 child
                                             # In this case there can be at most 3 children
                                             prob=c(0.5,0.5,0.5,0.5), 
                                             probConf=c(0,0,1), 
                                             probInts=c(0.5,0.5,0.5,0.5), 
                                             weights=0,
                                             blind=TRUE,
                                             IBlind=NULL,
                                             YBlind=NULL,
                                             require_path=FALSE,
                                             IFactor=1.0,
                                             shuffleObs=FALSE,
                                             debug=TRUE) {
  # Total number of variables (regime + confounders + intervention variables + observed )
  P <- 1 + numConf + numInts + p
  
  # If B is given, use that as adjacency matrix (p, numConf and numInts must still be set; no checking is done)
  # Otherwise, generate B
  if ( is.null(B) ) {
    B <- simul_interventions_NIPS.generateB(P, p, numConf, numInts, prob, probConf, probInts, weights, debug+FALSE)
  } else if ( !is.null(weights) ) {
    # If B is NULL and weights is not, leave B exactly as is.
    # Replace nonzeros in B by random weights (except for arrows R --> I_i)
    for (i in 2:P) {
      for (j in 1:P) {
        if ( B[i,j] ) {
          if( weights == 1 ) {
            randomweight <- runif(1, -sqrt(3), sqrt(3))
          } else if( weights == -1 ) {
            randomweight <- runif(1, 2./3., 4./3.) * (rbinom(Nchildren, 1, 0.5)*2-1)
          } else if( weights == -2 ) {
            randomweight <- runif(1, 0, 2*sqrt(3))
          } else {
            randomweight <- rnorm(1, 0.2,0.8) * sample(c(-1,1),1)
          }
          B[i,j] <- randomweight
        }
      }
    }
  }

  # (regime + confounders + intervention variables + observed ), P = total
  B[(numConf+2):(numConf+1+numInts), (numConf+numInts+2):P] <- IFactor * B[(numConf+2):(numConf+1+numInts), (numConf+numInts+2):P]
  
  # Apply random permutation of variables
  perm <- 1:p
  if ( shuffleObs ) {
    perm <- sample(1:p)
  }

  # Only observed variables are permuted
  permutedB <- B[(numConf+numInts+2):P,(numConf+numInts+2):P]
  permutedB <- permutedB[perm,perm]

  ### Convert to standard format data.
  
  gen_M <- list(G=array(0,c(p,p)), Ge=array(0,c(p,p)), Gs=array(0,c(p,p)))
  gen_M$perm <- perm
  gen_M$B <- t(permutedB)
  
  # No confounders
  fullperm <- c((numConf+numInts+1+perm), 1, (numConf+2):(numConf+numInts+1))
  gen_M$fullB <- t(B[fullperm,fullperm])
  gen_M$fullperm <- fullperm
  
  if ( blind && (is.null(IBlind) || is.null(YBlind)) ) {
    IBlind <- sample(x=numInts, size=1)
    YBlind <- sample(x=p, size=1)
    gen_M$IBlind <- IBlind
    gen_M$YBlind <- YBlind
    
    if (debug) cat("YBlind:", YBlind, "IBlind:", IBlind)
    
    # No direct effect
    gen_M$fullB[YBlind, p+1+IBlind] <- 0
  }
  
  
  gen_M$fullGe <- array(0,c(nrow(gen_M$fullB),nrow(gen_M$fullB)))
  
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
  
  if (blind && gen_M$fullG[gen_M$YBlind, p+1+gen_M$IBlind] == 1) {
    stop("Error, direct effect")
  }
  
  gen_M
}

simul_interventions_NIPS.generateModel <- function(B=NULL, p=4, numConf=2, numInts = 2,
                                             # Each value is the probability of the number of children, e.g. 0.5 is the prob of 1 child
                                             # In this case there can be at most 3 children
                                             prob=c(0.125,0.5,0.25,0.125),
                                             probConf=c(0,0,1),
                                             probInts=c(0.125,0.5,0.25,0.125),
                                             weights=0,
                                             blind=TRUE,
                                             IBlind=NULL,
                                             YBlind=NULL,
                                             filter_require_path_IY=FALSE,
                                             IFactor=1.0,
                                             shuffleObs=TRUE,
                                             debug=TRUE) {
  done <- FALSE
  while (!done) {
    all_filters_pass <- TRUE
    gen_M <- simul_interventions_NIPS.generateModel_unfiltered(B=B, p=p, numConf=numConf, numInts=numInts,
                                                               prob=prob, probConf=probConf, probInts=probInts, weights=weights,
                                                               blind=blind, IBlind=IBlind, YBlind=YBlind,
                                                               IFactor=IFactor, shuffleObs=shuffleObs, debug=debug)

    if (filter_require_path_IY) {
      # reject models in which there is no path from IBlind to YBlind
      tryIBlind <- gen_M$IBlind
      tryYBlind <- gen_M$YBlind
      if (gen_M$fullC[p+1+tryIBlind, tryYBlind] != 1) {
        all_filters_pass <- FALSE
        if (debug) {
          cat("Model rejected by filter_require_path_IY\n")
        }
      }
    }
    
    if (all_filters_pass) {
      done <- TRUE
    }
  }
  gen_M
}

simul_interventions_NIPS.generateData <- function(M, disturbance_distribution=function(x){0.1*rnorm(x,0,0.8)}, base_function=identity, 
                                            nObs=500, nPerInt=50, weights=0) {
  B <- M$original_B
  numConf <- M$numConf
  numInts <- M$numInts
  p <- nrow(M$G)

  IBlind <- M$IBlind
  YBlind <- M$YBlind
  
  P <- 1 + numConf + numInts + p

  nAll <- nObs + nPerInt*numInts

  ## Simulate observational data
  ##allObs <- matrix(rnorm(nObs*P),nrow=nObs) %*% diag(scales * noiseSigma)
  ##X <- allObs / (eye(P) - B)   # SLOW
  #obs <- t(forwardsolve(t(diag(1,P) - B), t(allObs), transpose=FALSE))
  
  # Simulate observational data
  alldata <- matrix(0, P, nrow=nAll)
  # Not very efficient: one by one.
  for (j in 1:nAll) {
    V <- array(0, c(P))
    # Choose regime randomly:
    #V[1] <- sample(x = numInts + 1, 1)
    # Choose regime according to nObs, nPerInt:
    if ( j <= nObs ) {
      V[1] <- 1
    } else {
      V[1] <- 1 + ceiling((j - nObs) / nPerInt)
    }
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
  
  # Apply random permutation of variables
  perm <- M$perm
  obs <- obs[,perm]

  regime <- alldata[,1]
  ints <- alldata[,(numConf+2):(numConf+numInts+1)]
  
  final_data <- data.frame(cbind(obs, regime))
  names(final_data) <- c(paste("v", 1:p, sep=""), "experiment")

  blinded_obs <- obs
  for (j in 1:nAll) {
    if ( !is.null(IBlind) && regime[j] == IBlind + 1 ) {
      blinded_obs[j,YBlind] <- NaN
    }
  }

  blinded_data <- data.frame(cbind(blinded_obs, regime))
  names(blinded_data) <- c(paste("v", 1:p, sep=""), "experiment")
  
  ### Convert to standard format data
  D <- list()
  #E <- experimentConfiguration(p, "passive")
  D[[1]] <- list()
  D[[1]]$M<-M
  #D[[1]]$e<-E[1,]
  D[[1]]$N<-nAll
  D[[1]]$all_data <- alldata # order of columns: regime, latent confounders, intervention vars, observed vars
  D[[1]]$data <- final_data # order of columns: observed vars, regime
  D[[1]]$data_blinded <- blinded_data
  D[[1]]$Cx<-cov(D[[1]]$data)
  
  
  list(M=M, D=D)
}
