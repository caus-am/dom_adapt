G2st <- function(G) {
  # Turns a graph G to a string representation, used to passed graphs to bnlearn.
  st<-c()
  for ( i in 1:nrow(G) ) {
    if ( any(G[i,]==1)) {
      st<-paste(st,'[X',i,'|',paste('X',which(G[i,]==1),collapse=':',sep=''),']',sep=''  )
    } else {
      st<-paste(st,'[X',i,']',sep='')
    }
  }
  st
}