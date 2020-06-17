dec.to.bin <- function( x, ndigits ) {
  # Changes a decimal integer x into a binary number with ndigits.

  Base.b <- array(NA, dim=c(1, ndigits))
  for(i in 1:ndigits){
    Base.b[, ndigits-i+1] <- (x %% 2)
    x <- (x %/% 2)
  }
  Base.b[1, ]
}


dec.to.bin.test<-function(ndigits=3) {
  # Test dec.to.bin  and bin.to.dec functions.
  for ( i in 0:(2^ndigits-1)) {
    cat(i,'->',dec.to.bin(i,ndigits),
	  '->',bin.to.dec(dec.to.bin(i,ndigits)),'\n')
  }
}

bin.to.dec <- function( Base.b ) {
  #Changes a binary number into a decimal integer.
  #REMEMBER TO +1 IF USED FOR INDEXING A BINARY CPT!

  ndigits = length(Base.b)
  sum(Base.b*2^((ndigits-1):0))
}



# # Changes a decimal integer x into a binary number with ndigits.
# 
# dec.to.bin  =  function( x, ndigits ) {
# 
#   Base.b  =  array(NA, dim=c(1, ndigits))
#   for(i in 1:ndigits){
#     Base.b[, ndigits-i+1]  =  (x %% 2)
#     x  =  (x %/% 2)
#   }
#   Base.b[1, ]
# }