createData<-function(M, e, samples=1000) {
  # Creates data.
  # INPUT:
  # M - model from generateModel.
  # e - experiment vector, e[i]=1 if x[i] intervened, 0 otherwise
  # samples - number of samples, possibly Inf.
  # OUTPUT:
  # D -list
  # D$data - sample data, NA if infinite sample data requested.
  # D$Cx   - Data covariance matrix.
  n <- nrow(M$B)

  D <- list()
  J <- diag(e) # interventions in this experiment.
  U <- diag(n)-J # I where 0 means an intervention.

  #first create the true covariance matrix
  # aka what would be observed in the infinite sample limit
  A<-mpinv(diag(n)-U%*%M$B)
  # mpinv matrix with 0 row where intervened.
  Cx<-A%*%(J+U%*%M$Ce%*%U)%*%t(A)

  # also store true model
  D$M<-M

  if ( is.infinite(samples) ) {
    D$data<-NA
    D$Cx<-Cx
  } else {   
    EIG<-eigen(Cx)
    sCx<-EIG$vectors%*%diag(sqrt(EIG$values))%*%t(EIG$vectors)
  
    X<-array(rnorm(samples*n),c(samples,n))
    D$data<-X%*%sCx

    D$Cx<-cov(D$data)
    #cat('Difference to true:',max(abs(D$Cx-Cx)),'\n')
    #browser()
  }
  D$e<-e
  D$N<-samples

  D
}
