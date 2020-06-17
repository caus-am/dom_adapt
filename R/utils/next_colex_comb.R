next_colex_comb<-function(x) {
#For a 0-1 vectors gives a next 
#Can be used to quickly iterate over all subset of a given size.
#Just start with (1,1,1,….,1,0,…,0). In the end returns 0.
#The underlying mechanism to determine the successor is to determine the lowest block of ones 
#and move its highest bit one position up. 
#The rest of the block is then moved to the low end of the word.
  j=1;
  for ( i in index(1,(length(x)-1)) ) {
    if ( x[i] == 1 ) {
      if ( x[i + 1] == 0 ) {
	x[i]<-0;x[i+1]=1;
	return(x); #switch bit to left
      } else {
	x[i]<-0;x[j]=1;j<-j+1;
      }
    }
  }
  return(NA)
}