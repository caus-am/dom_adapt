getm <- function(vars, C, n) {
  # Gives a integer representation of the M set given
  # the number of variables and conditioning set C.
  vec <- rep(1, n)
  #vec[x]<-vec[y]<-0
  vec[vars]<-0
  vec[C]<-0
  #bin.to.dec(vec)
  bin.to.dec(rev(vec))
}