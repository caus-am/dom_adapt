pattern2dag<-function(M) {
  #makes a dag out of a pdag by orienting unoriented edges
  #after this we can read of the independencies
  #first randomly orient any bidirected edges
  
  #make any bidirected edges as undirected so they will be directed
  M$Gs<-1*(M$Gs==1 | M$Ge==1)
  #notice that there is a problem
  
  #while ( any(M$Ge!= 0 ) ) {
  #  vars<-which(M$Ge == 1, arr.index=TRUE)[1,]
    #if ( sample(c("forward","backward"),1) == "forward" ) {
    #  M$G[vars[2],vars[1]]<-1
    #} else {
    #  M$G[vars[1],vars[2]]<-1      
    #}    
  #  M$Ge[vars[1],vars[2]]<-M$Ge[vars[2],vars[1]]<-0
  #}
  
  true_colliders<-unshielded_colliders( M$G, 1*(M$Gs==1 | M$G==1 | t(M$G==1)) )
  #print(true_colliders)
  #browser()
  #the remaining 
  require('combinat')
  for ( perm in permn(1:nrow(M$G)) ) {
    #cat('order:',perm,'\n')
    #other edges not ok with this order
    if ( any( (M$G[perm,perm] == 1) & upper.tri(M$G) ) ) next
    G<-M$G
    Gs<-M$Gs
    #browser()
    while ( any(Gs!= 0 ) ) {
      vars<-which(Gs == 1, arr.ind=TRUE)[1,]
      if ( which(perm == vars[1]) < which(perm == vars[2]) ) {
        G[vars[2],vars[1]]<-1
      } else {
        G[vars[1],vars[2]]<-1      
      }    
      Gs[vars[1],vars[2]]<-Gs[vars[2],vars[1]]<-0
    }
    #now should check the unshielded colliders
    colliders<-unshielded_colliders(G)
    if ( any ( dim(colliders) != dim(true_colliders)) ) next
    if ( any ( colliders != true_colliders ) ) next
    
    #if got this far everything OK
    return(list(G=G,Gs=array(0,dim(G)),Ge=array(0,dim(G))))
  }
  
  #stop("Could not find a DAG for the class.")
  
  NULL
}



