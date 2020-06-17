mpinv <- function(X)
{

	# Moore-Penrose inverse (by Thorsten Hothorn, downloaded from
	# http://tolstoy.newcastle.edu.au/R/help/99a/1118.htm)
	#
	# Example usage:
	# X <- cbind(1, diag(3)); # singular matrix
	# y <- 1:3
	# b2 <- mpinv(X)%*%y # least square fit using Moore-Penrose
	# X%*%b2 # == y 	
	#
	
	Eps <- 100 * .Machine$double.eps

	# singular value decomposition
	s <- svd(X)
	d <- s$d
	m <- length(d)
	if (!(is.vector(d)))
		return(t(s$v%*%(1/d)%*%t(s$u)))
	# remove eigenvalues equal zero
	d <- d[d > Eps]
	notnull <- length(d) 
	if (notnull == 1)
	{
		inv <- 1/d
	} else {
		#browser()
		#inv <- solve(diag(d))
		inv <- diag(1/d) #correction made by A. Hyttinen
	}
	# add rows, columns of zeros if needed 
	if (notnull != m)
	{
		inv <- cbind(inv, matrix(0, nrow=notnull, ncol=(m - notnull)))
		inv <- rbind(inv, matrix(0, nrow=(m-notnull), ncol=m))
	} 

	# compute Moore-Penrose
	mp <- s$v%*%inv%*%t(s$u)

	# set very small values to zero
	mp[abs(mp) < Eps] <- 0
	return(mp)
}

