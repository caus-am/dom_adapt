#A slightly modified version of the pcalg code.
#See http://cran.r-project.org/web/packages/pcalg/index.html for references.

pdsep_mod<-function (skel, suffStat, indepTest, p, sepset, pMax, NAdelete = TRUE, 
          verbose = FALSE, alpha, unfVect = NULL, biCC = FALSE) 
{
  #browser()
  G <- (as(skel, "matrix") != 0)
  n.edgetests <- rep(0, 1000)
  max.ord <- 0
  allPdsep <- vector("list", p)
  if (biCC) {
    conn.comp <- lapply(biConnComp(skel), as.numeric)
  }
  if (any(G)) {
    amat <- G
    amat[amat == TRUE] <- 1
    amat[amat == FALSE] <- 0
    ind <- which(amat == 1, arr.ind = TRUE)
    for (i in 1:dim(ind)[1]) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(amat[y, ] != 0), x)
      if (length(allZ) > 0) {
        for (j in seq_along(allZ)) {
          z <- allZ[j]
          if (amat[x, z] == 0 && !((y %in% sepset[[x]][[z]]) | 
                                     (y %in% sepset[[z]][[x]]))) {
            if (length(unfVect) == 0) {
              amat[x, y] <- amat[z, y] <- 2
              if (verbose) 
                cat(x, "o->", y, "<-o", z, "\n")
            }
            else {
              if (!any(unfVect == triple2numb(p, x, y, 
                                              z)) && !any(unfVect == triple2numb(p, 
                                                                                 z, y, x))) {
                amat[x, y] <- amat[z, y] <- 2
                if (verbose) 
                  cat(x, "o->", y, "<-o", z, "\n")
              }
            }
          }
        }
      }
    }
    for (x in 1:p) {
      if (any(x.has <- amat[x, ] != 0)) {
        allPdsep[[x]] <- qreach(x, amat)
        tf1 <- setdiff(allPdsep[[x]], x)
        if (verbose) {
          cat("Possible D-Sep of", x, " is: ", allPdsep[[x]], 
              "\n")
        }
        adj.x <- (1:p)[amat[x, ] != 0]
        for (y in adj.x) {
          tf <- setdiff(tf1, y)
          diff.set <- setdiff(tf, adj.x)
          if (biCC) {
            index.conncomp <- 0
            found.conncomp <- FALSE
            while ((!found.conncomp) & (index.conncomp < 
                                          length(conn.comp))) {
              index.conncomp <- index.conncomp + 1
              if (x %in% conn.comp[[index.conncomp]] && 
                    y %in% conn.comp[[index.conncomp]]) {
                found.conncomp <- TRUE
              }
            }
            bi.conn.comp <- setdiff(conn.comp[[index.conncomp]], 
                                    c(x, y))
            tmp.tf <- intersect(tf, bi.conn.comp)
            tf <- tmp.tf
            if (verbose) {
              cat("Possible D-Sep of", x, "and", y, "intersected with the biconnected component is", 
                  tf, "\n")
            }
          }
          if (length(diff.set) > 0) {
            done <- FALSE
            ord <- 0
            while (!done && ord < length(tf)) {
              ord <- ord + 1
              if (ord > max.ord) 
                max.ord <- ord
              if (ord == 1) {
                for (j in seq_along(diff.set)) {
                  pval <- indepTest(x, y, diff.set[j], 
                                    suffStat)
                  n.edgetests[ord + 1] <- n.edgetests[ord + 
                                                        1] + 1
                  if (is.na(pval)) 
                    pval <- if (NAdelete) 
                      1
                  else 0
                  if (pval > pMax[x, y]) 
                    pMax[x, y] <- pval
                  if (pval >= alpha) {
                    amat[x, y] <- amat[y, x] <- 0
                    sepset[[x]][[y]] <- sepset[[y]][[x]] <- diff.set[j]
                    done <- TRUE
                    if (verbose) 
                      cat("x=", x, " y=", y, " S=", diff.set[j], 
                          ": pval =", pval, "\n")
                    break
                  }
                }
              }
              else if (ord <= length(adj.x)) {
                tmp.combn <- combn(tf, ord)
                #Correction by antti h.
                if ( is.vector(tmp.combn)) tmp.combn<-array(tmp.combn,c(length(tmp.combn),1))
                #browser()
                for (k in 1:dim(tmp.combn)[2]) {
                  tmp.ii <- tmp.combn[, k] %in% adj.x
                  if (any(!tmp.ii)) {
                    pval <- indepTest(x, y, tmp.combn[, 
                                                      k], suffStat)
                    n.edgetests[ord + 1] <- n.edgetests[ord + 
                                                          1] + 1
                    if (is.na(pval)) 
                      pval <- if (NAdelete) 
                        1
                    else 0
                    if (pval > pMax[x, y]) 
                      pMax[x, y] <- pval
                    if (pval >= alpha) {
                      amat[x, y] <- amat[y, x] <- 0
                      sepset[[x]][[y]] <- sepset[[y]][[x]] <- tmp.combn[, 
                                                                        k]
                      done <- TRUE
                      if (verbose) {
                        cat("x=", x, " y=", y, " S=", 
                            tmp.combn[, k], ": pval =", 
                            pval, "\n")
                      }
                      break
                    }
                  }
                }
              }
              else {
                tmp.combn <- combn(tf, ord)
                #Correction by antti h.
                if ( is.vector(tmp.combn)) tmp.combn<-array(tmp.combn,c(length(tmp.combn),1))
                
                for (k in 1:dim(tmp.combn)[2]) {
                  pval <- indepTest(x, y, tmp.combn[, 
                                                    k], suffStat)
                  n.edgetests[ord + 1] <- n.edgetests[ord + 
                                                        1] + 1
                  if (is.na(pval)) 
                    pval <- if (NAdelete) 
                      1
                  else 0
                  if (pval > pMax[x, y]) 
                    pMax[x, y] <- pval
                  if (pval >= alpha) {
                    amat[x, y] <- amat[y, x] <- 0
                    sepset[[x]][[y]] <- sepset[[y]][[x]] <- tmp.combn[, 
                                                                      k]
                    done <- TRUE
                    if (verbose) 
                      cat("x=", x, " y=", y, " S=", tmp.combn[, 
                                                              k], ": pval =", pval, "\n")
                    break
                  }
                }
              }
            }
          }
        }
      }
    }
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
    G[amat == 2] <- TRUE
  }
  list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep, 
       max.ord = max.ord, n.edgetests = n.edgetests[1:(max.ord + 
                                                         1)])
}