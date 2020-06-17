directed_reachable <- function(x, y, C, J=c(), M, verbose=0) {
  # Implements a d-separation oracle.
  # x, y    variables considered
  # C      vector of variable indexes in the conditioning set
  # J      vector of variable indexes in the intervention set  
  # note M is model here, M$G, M$Ge, and M$Gs are all used.
  
  # Note that this is not currently the most efficient implementation.
  # -One could use matrix operations
  # -One could calculate more of the relations all at once.
  
  n <- nrow(M$G)
  
  # if asking for these there is some problem
  if ( x %in% C) browser()
  if ( y %in% C) browser()

  # for ( j in J ) { #take out the edges into intervened variables
  #  M$G[j,]<-0 
  #  M$Ge[j,]<-0
  #  M$Ge[,j]<-0
  # }
  M$G[J,] <- 0
  M$Ge[J,] <- 0
  M$Ge[,J] <- 0

  HH <- M$Ge #symmetric head head paths
  
  if ( is.null(M$Gs) ) { # there should generally be no head-head paths in the beginning
    TT <- array(0,c(n,n)) # this is for tail tail paths
  } else {
    TT <- M$Gs
  }
  TH <- M$G # notice here TH[x,y] = 1 iff path y->x
  
  # TH and HH self paths do not affect the d-connectedness so they can be ignored
  diag(TH) <- (-1)
  diag(HH) <- (-1)
  
  if (verbose) {
    cat('Beginning state:\n')
    cat("TH:\n")
    print(TH)
    cat("TT:\n")
    print(TT)
    cat("HH:\n")
    print(HH)  
  }
   
  for ( node in 1:n ) { # doing either conditioning or marginalizing to all variables
    if (node == x || node == y) 
      next; # skip variables that are in the final graph
    
    # gather up all different kinds of parents, or connected nodes
    thpa<-which( TH[node,] == 1)
    htpa<-which( TH[,node] == 1) #aka children
    hhpa<-which( HH[node,] == 1 | HH[,node] == 1)
    ttpa<-which( TT[node,] == 1)
    
    if ( !(node %in% C) ) {
      #the marginalization operation is more difficult
      if ( verbose ) cat('Marginalizing:',node,'\n')
      #i-->node-->j
      for ( i in thpa ) {
        for ( j in htpa ) {
          TH[j,i]<-1  
        }
      }
      
      #i--> node --- j
      for ( i in thpa ) {
        for ( j in ttpa ) {
          TT[i,j]<-TT[j,i]<-1  
        }
      }
      #i --> node <-- is not ok
      #i --> node <-> is not ok
      #####################################
      
      #i <-> node --> j
      for ( i in hhpa ) {
        for ( j in htpa ) {
          HH[j,i]<-HH[i,j]<-1  
        }
      }
      
      #i <-> node --- j
      for ( i in hhpa ) {
        for ( j in ttpa ) {
          TH[i,j]<-1  
        }
      }
      #i <-> node <-- is not ok
      #i <-> node <-> is not ok
      ######################################
      
      # i --- node --- j
      #tail tail parents connected
      for ( i in ttpa ) { #connects node to itself as well so tt-self cycle is inherited
        for ( j in ttpa ) {
          TT[i,j]<-TT[j,i]<-1  
        }
      }
      #i --- node --> j
      for ( i in ttpa ) { #connects node to itself as well so tt-self cycle is inherited
        for ( j in htpa ) {
          TH[j,i]<-1  
        }
      }    
      # i --- node <-> j done already
      # i --- node <-- j done already
      ##############################################
      
      # i <-- node --> j
      for ( i in htpa ) { #connects node to itself as well so tt-self cycle is inherited
        for ( j in htpa ) {
          HH[i,j]<-HH[j,i]<-1  
        }
      }    
      #i <-- node <-> j done already
      #for ( i in htpa ) { #connects node to itself as well so tt-self cycle is inherited
      #  for ( j in hhpa ) {
      #    HH[i,j]<-HH[j,i]<-1  
      #  }
      #}    
      # i --- node <-> j done already
      # i --- node <-- j done already
      
      
    } #if node not in C
    if ( node %in% C || TT[node,node] == 1 ) {
        #notice the simplicity here!
        #an unconditioned node with a selfloop actually allows through all traffic!!!
      if ( verbose ) cat('Conditioning:',node,'\n')
      #only three options
      # i--> node <--j
      for ( i in thpa ) {
        for ( j in thpa ) { #notice that this connects that parents to them selves as well
          TT[i,j]<-TT[j,i]<-1  
        }
      } 
      # i<-> node <->j
      #hh parents need to be connected by head head nodes
      for ( i in hhpa ) {
        for ( j in hhpa ) {
          HH[i,j]<-HH[j,i]<-1  
        }
      } 
      # i<-> node <--j
      #connecting hh parent to th parent
      for ( i in hhpa ) {
        for (j in thpa ) {
          TH[i,j]<-1
        }
      }
    } #if node in C or a tail to tail path
    #only tailtail cycles are relevant to d-connection
    diag(TH)<-(-1)
    diag(HH)<-(-1)
    
    #now take the node away
    TH[node,]<-TH[,node]<-HH[,node]<-HH[node,]<-TT[,node]<-TT[node,]<-(-1)    
    if ( verbose )  cat("TH:\n")
    if ( verbose )  print(TH)
    if ( verbose )  cat("TT:\n")
    if ( verbose )  print(TT)
    if ( verbose )  cat("HH:\n")
    if ( verbose ) print(HH) 
  }
  #the nodes are connected if any of the paths is present in the end, where
  #all variables 
  HH[x,y] ==1 | TH[x,y] ==1 |TH[y,x] == 1 | HH[y,x] == 1 | TT[y,x] == 1
}

directed_reachable_test<-function(case = 1) {
  #testing function for the oracle
  M<-list()
  if (case == 1 ) {#richardson graph 1
    n=4
    M$G<-array(0,c(n,n))  
    M$Ge<-array(0,c(n,n))
    M$G[3,1]<-M$G[4,2]<-M$G[3,4]<-M$G[4,3]<-1
    J<-c()
  } else   if (case == 2 ) {#richardson graph 1
    n=4
    M$G<-array(0,c(n,n))  
    M$Ge<-array(0,c(n,n))
    M$G[3,1]<-M$G[4,2]<-M$G[3,4]<-M$G[4,3]<-1
    M$Ge[3,4]<-M$Ge[4,3]<-1
    J<-c()
  } else if (case == 3) {
    n=4 
    M$G<-array(0,c(n,n))  
    M$Ge<-array(0,c(n,n))
    M$G[1,2]<-M$G[4,3]<-1
    M$Ge[4,2]<-1
    J<-c()
  }
  result <- ""
  
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i >= j ) next
      remaining<-setdiff(1:n,c(i,j))
      for ( csetindex in index(0,2^(n-2)-1) ) {
        C<-remaining[dec.to.bin(csetindex,length(remaining))==1]
        if (   directed_reachable(i,j,C,J,M,verbose=0) ) {
          result <- paste(result, i,'_N_',j,
                          ' |',paste(C,collapse=','),' ||',paste(J,collapse=','),'\n',sep='')
        } else {
          result <- paste(result, i,'_||_',j,
                          ' |',paste(C,collapse=','),' ||',paste(J,collapse=','),'\n',sep='')       
        }
        #browser()
      }      
    }
  }
  
  cat(result)
  
  if (case == 1){
    expected_result <- "1_||_2 | ||
    1_N_2 |4 ||
    1_N_2 |3 ||
    1_||_2 |3,4 ||
    1_N_3 | ||
    1_N_3 |4 ||
    1_N_3 |2 ||
    1_N_3 |2,4 ||
    1_N_4 | ||
    1_N_4 |3 ||
    1_N_4 |2 ||
    1_N_4 |2,3 ||
    2_N_3 | ||
    2_N_3 |4 ||
    2_N_3 |1 ||
    2_N_3 |1,4 ||
    2_N_4 | ||
    2_N_4 |3 ||
    2_N_4 |1 ||
    2_N_4 |1,3 ||
    3_N_4 | ||
    3_N_4 |2 ||
    3_N_4 |1 ||
    3_N_4 |1,2 ||
    "
  }
  if (case == 2){  
    expected_result <- "1_||_2 | ||
    1_N_2 |4 ||
    1_N_2 |3 ||
    1_N_2 |3,4 ||
    1_N_3 | ||
    1_N_3 |4 ||
    1_N_3 |2 ||
    1_N_3 |2,4 ||
    1_N_4 | ||
    1_N_4 |3 ||
    1_N_4 |2 ||
    1_N_4 |2,3 ||
    2_N_3 | ||
    2_N_3 |4 ||
    2_N_3 |1 ||
    2_N_3 |1,4 ||
    2_N_4 | ||
    2_N_4 |3 ||
    2_N_4 |1 ||
    2_N_4 |1,3 ||
    3_N_4 | ||
    3_N_4 |2 ||
    3_N_4 |1 ||
    3_N_4 |1,2 ||
    "
  }
  if (case == 3){    
    expected_result <- "1_N_2 | ||
    1_N_2 |4 ||
    1_N_2 |3 ||
    1_N_2 |3,4 ||
    1_||_3 | ||
    1_N_3 |4 ||
    1_||_3 |2 ||
    1_||_3 |2,4 ||
    1_N_4 | ||
    1_N_4 |3 ||
    1_||_4 |2 ||
    1_||_4 |2,3 ||
    2_||_3 | ||
    2_N_3 |4 ||
    2_||_3 |1 ||
    2_N_3 |1,4 ||
    2_N_4 | ||
    2_N_4 |3 ||
    2_N_4 |1 ||
    2_N_4 |1,3 ||
    3_N_4 | ||
    3_N_4 |2 ||
    3_N_4 |1 ||
    3_N_4 |1,2 ||
    "
  }

  if(result != expected_result) {
    cat("[Error]: result not as expected.")
  }
}
