checkTriple2 =  function (a, 
    b, c, nbrsA, nbrsC, sepsetA, sepsetC, suffStat, indepTest, 
    alpha, version.unf = c(NA, NA), maj.rule = FALSE, verbose = FALSE, m.max = Inf) 
{
  nr.indep <- 0
    stopifnot(length(version.unf) == 2, version.unf %in% 1:2)
    tmp <- if (version.unf[2] == 2) 
        (b %in% sepsetA || b %in% sepsetC)
    version <- 0
    if ((nn <- length(nbrsA)) > 0) {
        allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
        
        # Modified 
        row.has.too_many_ones <- apply(allComb, 1, function(x){length(which(x==1)) > m.max})
        allComb <- allComb[!row.has.too_many_ones, , drop=FALSE]
        
        for (i in 1:nrow(allComb)) {
            S <- nbrsA[which(allComb[i, ] != 0)]
            pval <- indepTest(a, c, S, suffStat)
            if (verbose) 
                cat("a: S =", S, " - pval =", pval, "\n")
            if (pval >= alpha) {
                nr.indep <- nr.indep + 1
                tmp <- c(tmp, b %in% S)
                version <- 1
            }
        }
    }
    if ((nn <- length(nbrsC)) > 0) {
        allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
        
        # Modified 
        row.has.too_many_ones <- apply(allComb, 1, function(x){length(which(x==1)) > m.max})
        allComb <- allComb[! row.has.too_many_ones, , drop=FALSE]
        
        for (i in 1:nrow(allComb)) {
            S <- nbrsC[which(allComb[i, ] != 0)]
            pval <- indepTest(a, c, S, suffStat)
            if (verbose) 
                cat("c: S =", S, " - pval =", pval, "\n")
            if (pval >= alpha) {
                nr.indep <- nr.indep + 1
                tmp <- c(tmp, b %in% S)
                version <- 1
            }
        }
    }
    if (version.unf[1] == 2 && nr.indep == 0) {
        version <- 2
    }
    if (is.null(tmp)) 
        tmp <- FALSE
    if (all(tmp)) {
        res <- 2
        if (!(b %in% sepsetA))
            sepsetA <- c(sepsetA, b)
        if (!(b %in% sepsetC))
            sepsetC <- c(sepsetC, b)
    }
    else {
        if (all(!tmp)) {
            res <- 1
            sepsetA <- setdiff(sepsetA, b)
            sepsetC <- setdiff(sepsetC, b)
        }
        else {
            if (!maj.rule) {
                res <- 3
            }
            else {
                if (sum(tmp)/length(tmp) < 0.5) {
                  res <- 1
                  sepsetA <- setdiff(sepsetA, b)
                  sepsetC <- setdiff(sepsetC, b)
                }
                else if (sum(tmp)/length(tmp) > 0.5) {
                  res <- 2
                  if (!(b %in% sepsetA))
                    sepsetA <- c(sepsetA, b)
                  if (!(b %in% sepsetC)) 
                    sepsetC <- c(sepsetC, b)
                }
                else if (sum(tmp)/length(tmp) == 0.5) {
                  res <- 3
                }
            }
        }
    }
    if (verbose && res == 3) 
        cat("Triple ambiguous\n")
    lapply(list(decision = res, version = version, SepsetA = sepsetA, 
        SepsetC = sepsetC), as.integer)
}