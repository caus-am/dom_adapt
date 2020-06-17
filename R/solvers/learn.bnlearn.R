learn.bnlearn<-function(D, n, p_threshold=1) {
  #Function for running the learning tasks for score-based learning.
  #This implements an exact search over the bnlearn package.
  #The function is rather a quick fix.
  #D    data set from createDataSet.
  #p_threshold - sparsity penalty for the BIC local score
  #
  #Results options (where to store the results and tmp files):
  #indFilename          - the serialized version of independences, default: "pipeline.ind"
  tic()
  df<-data.frame(D[[1]]$data)
  
  #exhaustive search, not optimally fast here since looping over orders
  #there does not seem to be a proper implementation of score-based learning
  #in the exact sense for R
  
  #empty graph
  #overall the best graph
  bestbestG<-array(0,c(n,n))
  bestbestscore<-(-Inf)
  ind<-0
  
  for ( corder in permn(n) ) { #loop through all causal orders
    ind<-ind+1
    #These are the best options for the order
    bestG<-array(0,c(n,n))
    bestscore<-bnlearn::score(as.bn( G2st(bestG) ),df,k=p_threshold*log(nrow(df))/2)
    for ( i in n:2 ) { #consider nodes from the end of the order to beginning
      G<-bestG #always start with the best here
      node<-corder[i]
        
      #cat('\tnode:',node,'\n')
        
      #find the possible parents for the node
      possible_parents<-corder[index(1,(i-1))]
        
      #print(possible_parents)
      #cat('\t\tpossible_parents:',possible_parents,'\n')
        
      bestscore<-(-Inf)
      for ( ipa in index(0,2^length(possible_parents)-1)) {#do not have consider 0
        G[node,possible_parents]<-dec.to.bin(ipa,length(possible_parents))
        score<-bnlearn::score(as.bn( G2st(G) ),df,k=p_threshold*log(nrow(df))/2)
        #cat(G[node,],'score=',score)
        if ( score > bestscore ) {
          bestG<-G
          bestscore<-score
          #cat('!')
        } 
        #cat('\n')
      }#parent configuration
    }#for i/node
    #now bestG is the best given this order, updating if it is the best given all orders
    if ( bestscore > bestbestscore) {
      bestbestG<-bestG
      bestbestscore<-bestscore
    }
    cat(ind,'bestscore:',bestscore,'bestbestscore:',bestbestscore,'\n')
  }#for order
  
  # make the L object from the best graph
  L <- list()
  L$G <- bestbestG
  L$Ge <- array(0, c(n,n))
  L$objective <- NA
  L$solving_time<-toc() 
  L$C <- dag2causes(L$G)
  
  L$G[which(L$G==0)] <- -1
  L$C[which(L$C==0)] <- -1
  
  # Note that these tested_independences are the ones from the oracle with schedule n-2.
  L
}


