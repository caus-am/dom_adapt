fci_mod2 = function (suffStat, 
    indepTest, alpha, labels, p, skel.method = c("stable", "original", 
        "stable.fast"), type = c("normal", "anytime", "adaptive"), 
    fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, 
    pdsep.max = Inf, rules = rep(TRUE, 10), doPdsep = TRUE, biCC = FALSE, 
    conservative = FALSE, maj.rule = FALSE, verbose = FALSE) 
{
    cl <- match.call()
    if (!missing(p)) 
        stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
            1, p >= 2)
    if (missing(labels)) {
        if (missing(p)) 
            stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
    }
    else {
        stopifnot(is.character(labels))
        if (missing(p)) {
            p <- length(labels)
        }
        else if (p != length(labels)) 
            stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
        else message("No need to specify 'p', when 'labels' is given")
    }
    type <- match.arg(type)
    if (type == "anytime" && m.max == Inf) 
        stop("To use the Anytime FCI you must specify a finite 'm.max'.")
    if (type == "adaptive" && m.max != Inf) 
        stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")
    if (conservative && maj.rule) 
        stop("Choose either conservative FCI or majority rule FCI")
    cl <- match.call()
    if (verbose) 
        cat("Compute Skeleton\n================\n")
    skel <- pcalg::skeleton(suffStat, indepTest, alpha, labels = labels, 
        method = "stable", fixedGaps = fixedGaps, fixedEdges = fixedEdges, 
        NAdelete = NAdelete, m.max = m.max, verbose = verbose)
    skel@call <- cl
    G <- as(skel@graph, "matrix")
    sepset <- skel@sepset
    pMax <- skel@pMax
    n.edgetestsSKEL <- skel@n.edgetests
    max.ordSKEL <- skel@max.ord
    allPdsep <- NA
    tripleList <- NULL
    if (doPdsep) {
        if (verbose) 
            cat("\nCompute PDSEP\n=============\n")
        pc.ci <- pc.cons.intern2(skel, suffStat, indepTest, alpha = alpha, 
            version.unf = c(1, 1), verbose = verbose,  maj.rule = FALSE, m.max= m.max)

        pdsepRes <- pcalg::pdsep(skel@graph, suffStat, indepTest = indepTest, 
            p = p, sepset = pc.ci$sk@sepset, alpha = alpha, pMax = pMax, 
            m.max = if (type == "adaptive") 
                max.ordSKEL
            else m.max, pdsep.max = pdsep.max, NAdelete = NAdelete, 
            unfVect = pc.ci$unfTripl, biCC = biCC, verbose = verbose)
        G <- pdsepRes$G
        sepset <- pdsepRes$sepset
        pMax <- pdsepRes$pMax
        allPdsep <- pdsepRes$allPdsep
        n.edgetestsPD <- pdsepRes$n.edgetests
        max.ordPD <- pdsepRes$max.ord
        if (conservative || maj.rule) {
            if (verbose) 
                cat("\nCheck v-structures conservatively\n=================================\n")
            tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), 
                call = cl, n = integer(0), max.ord = as.integer(max.ordSKEL), 
                n.edgetests = n.edgetestsSKEL, sepset = sepset, 
                pMax = pMax, zMin = matrix(NA, 1, 1))
            sk. <- pc.cons.intern2(tmp.pdsep, suffStat, indepTest, 
                alpha, verbose = verbose, version.unf = c(1, 
                  1), maj.rule = maj.rule, m.max= m.max)
            tripleList <- sk.$unfTripl
            sepset <- sk.$sk@sepset
            if (length(tripleList) >= 1 ){
              cat("\nUnfaithful triple\n")
            }
        }
    }
    else {
        n.edgetestsPD <- 0
        max.ordPD <- 0
        allPdsep <- vector("list", p)
        if (conservative || maj.rule) {
            if (verbose) 
                cat("\nCheck v-structures conservatively\n=================================\n")
            nopdsep <- pc.cons.intern(skel, suffStat, indepTest, 
                alpha, verbose = verbose, version.unf = c(2, 
                  1), maj.rule = maj.rule)
            tripleList <- nopdsep$unfTripl
            sepset <- nopdsep$sk@sepset
            if (length(tripleList) >= 1 ){
              cat("\nUnfaithful triple\n")
            }
        }
    }
    if (verbose) 
        cat("\nDirect egdes:\n=============\nUsing rules:", which(rules), 
            "\nCompute collider:\n")
    res <- udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList, 
        verbose = verbose)
    colnames(res) <- rownames(res) <- labels
    new("fciAlgo", amat = res, call = cl, n = integer(0), max.ord = as.integer(max.ordSKEL), 
        max.ordPDSEP = as.integer(max.ordPD), n.edgetests = n.edgetestsSKEL, 
        n.edgetestsPDSEP = n.edgetestsPD, sepset = sepset, pMax = pMax, 
        allPdsep = allPdsep)
}