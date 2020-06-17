build_tree <- function(aspSetsFullPath, n, tested_independences, verbose=0) {
  # Builds the encoding DAG and save it to a filename.
  require('hash')
  H <- hash()
  aspSetsFile <- file(aspSetsFullPath, "a")
  
  H[["0,0,0"]] <- 1 #this has been connected downward already.
  
  # Function that connects a key downward
  connect<-function(key) {
    if ( !is.null(H[[key]]) ) { #key already connected downward.
      if ( verbose ) cat('%(', key, '). STOP.\n', file = aspSetsFile)
      return ( 0 );
    }
    if ( verbose ) cat('%(', key, ').\n', file = aspSetsFile)
    conf<-as.numeric( strsplit(key,',')[[1]] )
    C<-which( rev(dec.to.bin(conf[1],n)) == 1 )
    J<-which( rev(dec.to.bin(conf[2],n)) == 1 )
    M<-which( rev(dec.to.bin(conf[3],n)) == 1 )

    #now connecting the J,C,M combination downward.
    if ( length(C) >= 1) { #first take out the conditioning
      set<-C
    } else if (length(M) >= 1 ) { #then get rid of marginalization
      set<-M
    } else { #finally take out intervention
      set<-J
    }
    
    #set<-union(union(C,M),J)
    if ( length(set) == 1 ) {
      z<-set[1]
    } else {
      z<-sample(set,1) #take a 
    }
    
    subconf<-conf
    if ( z %in% C ) { #notice here that if z is in C and J, it must be unconditioned first
      subconf[1]<-subconf[1]-2^(z-1)      
      cat('condition(',subconf[1],',',z,','
          ,conf[1],',',conf[2],',',conf[3],').\n',sep='', file = aspSetsFile)
      #cat("c",subconf[1],'j',conf[2],"m",conf[3],"->",
      #    "c",conf[1],'j',conf[2],"m",conf[3],";\n",sep="")
    } else if ( z %in% M ) {
      subconf[3]<-subconf[3]-2^(z-1)      
      cat('marginalize(',conf[1],',',conf[2],',',
          subconf[3],',',z,',',conf[3],').\n',sep='', file = aspSetsFile)   
      #cat("c",conf[1],'j',conf[2],"m",subconf[3],"->",
      #    "c",conf[1],'j',conf[2],"m",conf[3],";\n",sep="")
    } else if ( z %in% J ) {
      subconf[2]<-subconf[2]-2^(z-1)      
      cat('intervene(',conf[1],',',subconf[2],',',z,',',conf[2],
          ',',conf[3],').\n',sep='', file = aspSetsFile)
      #cat("c",conf[1],'j',subconf[2],"m",conf[3],"->",
      #    "c",conf[1],'j',conf[2],"m",conf[3],";\n",sep="")
    }
    #connect this to a lower one
    #cat( '%to',subconf,'.\n')
    H[[key]]<-1 #the key is connected
    subkey<-paste(subconf,collapse=',')
    #browser()
    #call connect to the subkey
    connect(subkey);
  }
  
  # Make sure all indeps are connected downward.
  i <- 0
  for ( indep in tested_independences) { #all indeps talk about different
    i <- i + 1
    if ( verbose ) 
      cat('%independence_index(', i, ').\n', file = aspSetsFile)
    connect( paste(indep$cset,indep$jset,indep$mset,sep=',') )
  }
  close(aspSetsFile)
  clear(H)
  invisible(0)
}

build_full<-function(aspSetsFullPath, n, tested_independences, verbose=TRUE) {
#Builds a full encoding DAG:
  #must build a tree from 0,0,0 so that all 
  
  require('hash')
  H<-hash()
  H[["0,0,0"]]<-1 #this has been connected
  
  aspSetsFile <- file(aspSetsFullPath, "a")
  
  connect<-function(key) {
    if ( !is.null(H[[key]]) ) {
      if ( verbose ) cat('%(', key, '). STOP.\n', file=aspSetsFile)
      return ( 0 );
    }
    if ( verbose ) cat('%(', key, ').\n', file=aspSetsFile)
    conf<-as.numeric( strsplit(key,',')[[1]] )
    C<-which( rev(dec.to.bin(conf[1],n)) == 1 )
    J<-which( rev(dec.to.bin(conf[2],n)) == 1 )
    M<-which( rev(dec.to.bin(conf[3],n)) == 1 )
    set<-union(union(C,M),J)
    for ( z in set ) {
      subconf<-conf
      if ( z %in% C ) {
        subconf[1]<-subconf[1]-2^(z-1)      
        cat('condition(',subconf[1],',',z,','
          ,conf[1],',',conf[2],',',conf[3],').\n',sep='', file=aspSetsFile)
      } else if ( z %in% M ) {
        subconf[3]<-subconf[3]-2^(z-1)      
        cat('marginalize(',conf[1],',',conf[2],',',
            subconf[3],',',z,',',conf[3],').\n',sep='', file=aspSetsFile)   
      } else if ( z %in% J ) { #if z is both in J and C, one must first do uncondition or unmarginalize
        subconf[2]<-subconf[2]-2^(z-1)      
        cat('intervene(',conf[1],',',subconf[2],',',z,',',conf[2],
          ',',conf[3],').\n',sep='', file=aspSetsFile)
      }
      #connect this to a lower one
      #cat( '%to',subconf,'.\n')
      H[[key]]<-1 #the key is connected
      subkey<-paste(subconf,collapse=',')
      #browser()
      #call connect to the subkey
      connect(subkey);
    }#for z in set
  }
  i<-0
  for ( indep in tested_independences) { #all indeps talk about different
    i<-i+1
    if ( verbose ) cat('%independence_index(',i,').\n',file=aspSetsFile)
    connect( paste(indep$cset,indep$jset,indep$mset,sep=',') )
  }
  
  close(aspSetsFile)
  clear(H)
  invisible(0)
}


############# DEAD CODE ######################################
# Keeping in case it helps understanding the one above.
#
build_tree.old<-function(aspSetsFullPath, n, tested_independences, verbose=0) {
  #must build a tree from 0,0,0 so that all 
  
  require('hash')
  H<-hash()
  
  H[["0,0,0"]]<-1 #this has been connected
  
  #sink('tmp/tree.pl')
  sink(aspSetsFullPath,append=TRUE)
  
  
  connect<-function(key) {
    if ( !is.null(H[[key]]) ) {
      if ( verbose ) cat('%(', key, '). STOP.\n')
      return ( 0 );
    }
    if ( verbose ) cat('%(', key, ').\n')
    conf<-as.numeric( strsplit(key,',')[[1]] )
    C<-which( rev(dec.to.bin(conf[1],n)) == 1 )
    J<-which( rev(dec.to.bin(conf[2],n)) == 1 )
    M<-which( rev(dec.to.bin(conf[3],n)) == 1 )
    set<-union(union(C,M),J)
    if ( length(set) == 1 ) {
      z<-set[1]
    } else {
      z<-sample(c(C,M,J),1) #take a 
    }
    subconf<-conf
    if ( z %in% C ) {
      subconf[1]<-subconf[1]-2^(z-1)      
      cat('condition(',subconf[1],',',z,','
          ,conf[1],',',conf[2],',',conf[3],').\n',sep='')
    } else if ( z %in% M ) {
      subconf[3]<-subconf[3]-2^(z-1)      
      cat('marginalize(',conf[1],',',conf[2],',',
          subconf[3],',',z,',',conf[3],').\n',sep='')      
    } else if ( z %in% J ) {
      subconf[2]<-subconf[2]-2^(z-1)      
      cat('intervene(',conf[1],',',subconf[2],',',z,',',conf[2],
          ',',conf[3],').\n',sep='')
    }
    #connect this to a lower one
    #cat( '%to',subconf,'.\n')
    H[[key]]<-1 #the key is connected
    subkey<-paste(subconf,collapse=',')
    #browser()
    #call connect to the subkey
    connect(subkey);
  }
  i<-0
  for ( indep in tested_independences) { #all indeps talk about different
    i<-i+1
    if ( verbose ) cat('%independence_index(',i,').\n')
    connect( paste(indep$cset,indep$jset,indep$mset,sep=',') )
  }
  
  sink()
  clear(H)
  invisible(0)
}