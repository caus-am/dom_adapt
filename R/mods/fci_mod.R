#Slightly modified version from the package pcalg, returns NULL if ambiguous triples found.
#See http://cran.r-project.org/web/packages/pcalg/index.html for references.

fci_mod<-function (suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL, 
          fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, rules = rep(TRUE, 10), doPdsep = TRUE, conservative = c(FALSE, FALSE), 
          biCC = FALSE, cons.rules = FALSE, labels = NA) 
{
  if (all(!is.na(labels))) {
    stopifnot(length(labels) == p)
  }
  else {
    labels <- as.character(1:p)
  }
  cl <- match.call()
  if (verbose) 
    cat("Compute Skeleton\n================\n")
  skel <- pcalg::skeleton(suffStat=suffStat, indepTest=indepTest, p=p, alpha=alpha, verbose = verbose, 
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges, NAdelete = NAdelete, 
                   m.max = m.max)
  G <- (as(skel@graph, "matrix") != 0)
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL
  if (conservative[1]) {
    sk <- pc.cons.intern(skel, suffStat, indepTest, alpha, 
                         verbose = verbose, version.unf = c(1, 2))
    tripleList <- sk$unfTripl
    #added 
    if ( conservative[1] && length(tripleList) >= 1 ) return(NULL)
  }
  if (doPdsep) {
    if (verbose) {
      cat("\nCompute PDSEP\n=============\ncompute collider...done\n")
    }
    pdsepRes <- pdsep_mod(skel@graph, suffStat, indepTest, p, 
                      sepset, pMax, NAdelete, verbose = verbose, alpha, 
                      unfVect = tripleList, biCC = biCC)
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    if (conservative[2]) {
      colnames(G) <- rownames(G) <- labels
      Gobject <- as(G, "graphNEL")
      tmp.pdsep <- new("pcAlgo", graph = Gobject, call = cl, 
                       n = integer(0), max.ord = as.integer(max.ordSKEL), 
                       n.edgetests = n.edgetestsSKEL, sepset = sepset, 
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk.pdsep <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, 
                                 alpha, verbose = verbose, version.unf = c(1, 
                                                                           2))
      tripleList <- sk.pdsep$unfTripl
      #added 
      if ( conservative[1] && length(tripleList) >= 1 ) return(NULL)
      
    }
  }
  else {
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
  }
  if (sum(G) == 0) {
    Gobject <- new("graphNEL", nodes = labels)
  }
  else {
    colnames(G) <- rownames(G) <- labels
    Gobject <- as(G, "graphNEL")
  }
  tmp <- new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
             max.ord = as.integer(max.ordSKEL), n.edgetests = n.edgetestsSKEL, 
             sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
  if (verbose) 
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules), 
        "\nCompute collider:\n")
  
  
  #added 
  if ( conservative[1] && length(tripleList) > 1 ) return(NULL)
  
  
  res <- if (numEdges(tmp@graph) > 0) 
    udag2pag(gInput = tmp, rules = rules, verbose = verbose, 
             unfVect = if (cons.rules) 
               tripleList)
  else G
  
  
  new("fciAlgo", amat = res, call = cl, n = integer(0), max.ord = as.integer(max.ordSKEL), 
      max.ordPDSEP = as.integer(max.ordPD), n.edgetests = n.edgetestsSKEL, 
      n.edgetestsPDSEP = n.edgetestsPD, sepset = sepset, pMax = pMax, 
      allPdsep = allPdsep)
}
