############# DEAD CODE ######################################
# Not used anywhere, keeping in case it helps debugging other code.
#
#functions for converting deal-like strings to a B matrix and back
dealstring2B<-function(dealstring,nvars) {
  varnames<-paste("X",c(1:nvars),sep="")#the names of the variables in order!
  B<-array(0,c(nvars,nvars))
  s<-substr(dealstring,start=2,stop=nchar(dealstring)-1) 
  #first rip off the beginning [ and ending ]
  s<-strsplit(s,"\\]\\[")[[1]]#the split with ][

  for (i in 1:length(s)) {
    si<-s[i]
    k<-strsplit(si,"\\|")[[1]]
    node<-which(varnames==k[1])
    if ( length(k) > 1 ) {
      parentsstring<-k[2]
    } else {
      parentsstring<-""
    }
    parents<-c()
    while( nchar(parentsstring) != 0 ) {
      for (j in 1:length(varnames)) {
        if ( substr(parentsstring,start=1,
                    stop=nchar(varnames[j])) == varnames[j] ) {
          parents<-c(parents,j)
          parentsstring<-substr(parentsstring,start=nchar(varnames[j])+2,
                                stop=nchar(parentsstring))
        }
      }
    }
    B[node,parents]<-1
  }
  B
}

B2dealstring<-function(B) {
  st<-''
  for ( i in 1:nrow(B) ) {
    pa<-which(B[i,]==1)
    if ( length(pa) == 0 ) {
      st<-paste(st,'[X',i,']',sep='')
    } else {
      st<-paste(st,'[X',i,'|',paste('X',pa,collapse=':',sep=''),']',sep='')
    }
  }
  st
}	
