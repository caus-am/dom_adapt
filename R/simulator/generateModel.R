generateModel <- function(n = 6, restrict = c('acyclic'), 
                        pedge = 1/(n-1), topology='random',
                        model = NULL,
                        samples = NULL,
                        confounder_proportion = 0.5,
                        verbose = 1) {
  # Function for generating random graph/linear gaussian model.
  # n  - number of observed variables
  # restrict - a vector with 'sufficient','acyclic' if these are wanted from the model.
  # pedge - the probability of an edge, this is good to be a decreasing function of n
  #         such that the density does not explode with increasing n
  #         (also note that for acyclic models the code , such that it is similar to )
  if (!is.null(samples) & is.null(model) ) {
    # No model necessary.
    return (NULL)
  }
  
  if (!is.null(model) ) {
    # No need to create a new model.
    if (verbose) {
      cat(" - Model already defined, no need to generate a new one.\n")
    }
    # Add causal relations to the model.
    model$C <- dag2causes(model$G)
    return (model)
  } else {
    if (verbose) {
      cat(" - Generating new model.", "\n")
    }
    # If the graph does not need to be sufficient.
    if (!any(restrict == 'sufficient') && topology != 'example2') {
      # Create a sufficient model for the causal part using the previous code.
      gen_M <- generateModel(n, union(restrict,"sufficient"), topology=topology, pedge=pedge*(1-confounder_proportion), verbose=verbose)
      
      gen_M$trueCe <- gen_M$Ce
      
      # Draw a sufficient model for the covariances.
      gen_Me <- generateModel(n, union(restrict,"sufficient"), topology=topology, pedge=pedge*confounder_proportion, verbose=verbose)
  
      # Then replace the covariance matrix A with the passively observed cov. matrix
      # of the error model: pseudo-inverse of ( I - B of the covariances)
      A <- mpinv(diag(n) - gen_Me$B)
      
      # Adding here some noise s.t. e are not deterministic functions
      # for example if we create noises e_x->e_y->e_z and the graph
      # does not have parents, x,y,z are is a deterministic function of es, 
      # and we will have x _||_z | y!
      # Ce_suff = A * Ce_cov * A' + diagonal n samples of N (0.5, 0.1)
      gen_M$Ce <- A%*%gen_Me$Ce%*%t(A)+diag( abs(0.5+0.1*rnorm(n)) )
      
      # Remember to update gen_M$Ge (confounders adjaciency matrix)
      gen_M$Ge <- abs(gen_M$Ce) > 1e-3    #replace also Ge which includes the true graph
      gen_M$Ce[abs(gen_M$Ce) < 1e-10] <- 0
      # Fix the diagonal to 0.
      diag(gen_M$Ge) <- 0
      return(gen_M)
    }
    
    # If the graph needs to be sufficient, sample the model.
    # Initialize the edges (G), confounders (Ge) and selection bias (Gs) adjanciency matrices.
    gen_M <- list(G=array(0,c(n,n)), Ge=array(0,c(n,n)), Gs=array(0,c(n,n)))
    
    if ( topology=='random') {
      #sample G either acyclic or possibly cyclic
      if ( any( restrict == 'acyclic') ) {
        order<-sample(1:n)
        pedge<-pedge*2 #multiply by two since fewer options
        gen_M$G<-array(sample(c(0,1),n*n,replace=TRUE,prob=c(1-pedge,pedge)),
                    c(n,n))*lower.tri(gen_M$G)
        diag(gen_M$G)<-0    
        order<-sample(1:n)
        gen_M$G[order,order]<-gen_M$G #sample the order again
      } else {
        # The generated model G is not acyclic.
        gen_M$G<-array( sample(c(0,1),n*n,replace=TRUE,prob=c(1-pedge,pedge)), c(n,n))
        diag(gen_M$G)<-0      
      }
      gen_M$B<-gen_M$G*matrix(runif(n*n,0.2,0.8),n,n)*matrix(sample(c(-1,1),n*n,replace=TRUE),n,n)  
  
    } else if (topology == 'example') {
      # The example is the y-structure.
      gen_M$G[3,1] <- gen_M$G[3,2] <- gen_M$G[4,3] <-1
      
      gen_M$B<- gen_M$G * matrix(runif(n*n,0.2,0.8),n,n) * 
        matrix(sample(c(-1,1),n*n,replace=TRUE),n,n)
    } else if (topology == 'example2') {
      #n <- 5
      n <- 4
      # 3 -> 1 and 2
      # 5 -> 2
      # 1 -> 4
      # 2 -> 4 (all inverted)
      gen_M$G[1,3] <- gen_M$G[2,3] <- gen_M$G[4,1] <- gen_M$G[4,2] <- gen_M$G[2,5] <-1
      
      gen_M$B<- gen_M$G * matrix(runif(n*n,0.2,0.8),n,n) * 
        matrix(sample(c(-1,1),n*n,replace=TRUE),n,n)
    }
  
    # Sample a diagonal for covariance matrices, the insufficiency part is handled earlier.
    # Sample n samples: |normal with mean 1 and variance 0.1|
    gen_M$Ce <- diag( abs(1 + 0.1*rnorm(n)) )
    
    # Add causal relations to the model.
    gen_M$C <- dag2causes(gen_M$G)
  
    if (verbose > 10) {
      cat("True G:\n")
      print(gen_M$G)
      cat("True Ge:\n")
      print(gen_M$Ge)
      cat("True C (transitive closure of G):\n")
      print(gen_M$C)
      cat("Coefficients of G:\n")
      print(gen_M$B)
      cat("Covariance of G:\n")
      print(gen_M$Ce)
    }
    gen_M
  }
}
