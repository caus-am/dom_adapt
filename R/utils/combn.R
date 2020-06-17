combn<-function (x, m, fun = NULL, simplify = TRUE, ...) 
{
  res<-combinat::combn(x,m,fun,simplify,...) 
  #this corrects the annoying R functionality when 
  #want 2 combinations of 2 elements, combn returns a vector
  if (is.vector(res)) res<-array(res,c(length(res),1))
  #browser()
  res
}