unshielded_colliders<-function( G, skel=1*(G==1 | t(G==1)) ) {
  #function that finds unshielded colliders
  #skel variable is used to check the shields
  #Note that the order of unshielded colliders is constant,
  #so it is enough the tests the colliders array to see whether
  #two graphs have the same set of colliders.
  colliders<-array(0,c(0,3))
  for ( i in 1:nrow(G) ) { #the middle node
    for ( j in 1:nrow(G) ) {
      if ( j == i ) next
      for ( k in 1:nrow(G) ){
        if ( k == i ) next
        if ( k <= j ) next #k==j we do not want, also if k <= j...
        #if j->i<-k and no shield 
        if ( G[i,j]==1 && G[i,k]==1 && skel[j,k] == 0 ) {
          colliders<-rbind(colliders,c(j,i,k))
        }
      }
    }
  }
  colliders
}
