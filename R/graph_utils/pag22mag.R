pag22mag<-function(M) {
  #Turns a pag into a mag, such that we can read off the independencies from that.
  #One simply cannot read off indeps from a pag, on has to first turn it into a mag,
  #and then read of the indeps.
  
  #First look at the true_colliders.
  true_colliders<-unshielded_colliders( 1*(M$G | M$Ge),#incoming arcs, and skeleton to show the shields 
                                        1*(M$Gs==1 | M$Ge==1 | M$G==1 | t(M$G==1) | M$Gcircles==1 ) )

  for ( perm in permn(1:nrow(M$G)) ) {#go through all permutations
    
    #if edges not ok with this order, try next
    if ( any( (M$G[perm,perm] == 1) & upper.tri(M$G) ) ) next
    
    L<-M
    
    #direct all circles according to the order
    while ( any(L$Gcircles!= 0 ) ) {
      vars<-which(L$Gcircles == 1, arr.ind=TRUE)[1,]
      if ( which(perm == vars[1]) < which(perm == vars[2]) ) {
        L$G[vars[2],vars[1]]<-1
      } else {
        L$G[vars[1],vars[2]]<-1      
      }    
      L$Gcircles[vars[1],vars[2]]<-L$Gcircles[vars[2],vars[1]]<-0
    }
    
    #get the colliders
    colliders<-unshielded_colliders( 1*(L$G | L$Ge), 
                                     1*(L$Gs==1 | L$Ge==1 | L$G==1 | t(L$G==1) ) )

    #if there are no new unshielded colliders, the orientation was successful
    if ( any ( dim(colliders) != dim(true_colliders)) ) next
    if ( any ( colliders != true_colliders ) ) next
    
    #if got this far everything OK
    L$Gcircles<-NULL
    return(L)
  }
  
  #stop("Could not find a DAG for the class.")
  
  NULL
}

