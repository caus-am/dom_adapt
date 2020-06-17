#A slightly modified version of the pcalg code.
#See http://cran.r-project.org/web/packages/pcalg/index.html for references.

pc.cons.intern_mod<-function (sk, suffStat, indepTest, alpha, verbose = FALSE, version.unf = c(NA, NA)) 
{
  skelObj <- list(G = as(sk@graph, "matrix"), sepset = sk@sepset)
  g <- skelObj$G
  stopifnot(all(g == t(g)))
  p <- dim(g)[1]
  unfTripl <- vers <- rep(NA, min(p * p, 1e+05))
  #added here by Antti 
  counter<-0
  if (sum(skelObj$G) > 0) {
    ind <- which(g == 1, arr.ind = TRUE)
    tripleMatrix <- NULL
    for (i in 1:dim(ind)[1]) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      allC <- setdiff(which(g[b, ] == 1), a)
      newC <- allC[g[a, allC] == 0]
      tmpMatrix <- cbind(rep(a, length(newC)), rep(b, length(newC)), 
                         newC)
      tripleMatrix <- rbind(tripleMatrix, tmpMatrix)
      colnames(tripleMatrix) <- c("", "", "")
    }
    #browser()
    if (dim(tripleMatrix)[1] >= 1) {
      deleteDupl <- rep(0, dim(tripleMatrix)[1])
      for (i in 1:dim(tripleMatrix)[1]) {
        if (tripleMatrix[i, 1] > tripleMatrix[i, 3]) {
          deleteDupl[i] <- 1
        }
      }
      deletedupl <- which(deleteDupl == 1)
      if (length(deletedupl) > 0) {
        finalMatrix <- tripleMatrix[-deletedupl, , drop = FALSE]
      }
      else {
        finalMatrix <- tripleMatrix
      }
      counter <- 0
      #browser()
      if (dim(finalMatrix)[1] >= 1) {
        for (i in 1:dim(finalMatrix)[1]) {
          if (counter == (length(unfTripl) - 1)) {
            tmp.vec <- rep(NA, min(p * p, 1e+05))
            unfTripl <- c(unfTripl, tmp.vec)
            vers <- c(vers, tmp.vec)
          }
          a <- finalMatrix[i, 1]
          b <- finalMatrix[i, 2]
          c <- finalMatrix[i, 3]
          nbrsA <- which(g[, a] != 0)
          nbrsC <- which(g[, c] != 0)
          if (verbose) {
            cat("Sepset by skelet:", skelObj$sepset[[a]][[c]], 
                "and", skelObj$sepset[[c]][[a]], "\n")
          }
          resTriple <- checkTriple(a, b, c, nbrsA, nbrsC, 
                                   skelObj$sepset[[a]][[c]], skelObj$sepset[[c]][[a]], 
                                   suffStat, indepTest, alpha, verbose = verbose, 
                                   version.unf = version.unf)
          if (resTriple$decision == 3) {
            counter <- counter + 1
            unfTripl[counter] <- triple2numb(p, a, b, 
                                             c)
            vers[counter] <- resTriple$version
            if (verbose) {
              cat("new unfTriple:", a, b, c, "\n")
            }
          }
          if ((version.unf[1] == 2) && (resTriple$version == 
                                          2) && (resTriple$decision != 3)) {
            counter <- counter + 1
            unfTripl[counter] <- triple2numb(p, a, b, 
                                             c)
            vers[counter] <- resTriple$version
            if (verbose) {
              cat("new unfTriple:", a, b, c, "\n")
            }
          }
        }
      }
    }
  }
  #browser()
  length(unfTripl) <- length(vers) <- counter
  list(unfTripl = unfTripl, vers = vers)
}