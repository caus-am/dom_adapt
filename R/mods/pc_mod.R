#Slightly modified version from the package pcalg, returns NULL if ambiguous triples found
#See http://cran.r-project.org/web/packages/pcalg/index.html for references.

pc_mod<-function (suffStat, indepTest, p, alpha, verbose = FALSE, fixedGaps = NULL, 
          fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, u2pd = "rand", conservative = FALSE) 
{
  cl <- match.call()
  skel <- pcalg::skeleton(suffStat=suffStat, indepTest=indepTest, p=p, alpha=alpha, verbose = verbose, 
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges, NAdelete = NAdelete, 
                   m.max = m.max)
  if (!conservative) {
    switch(u2pd, rand = udag2pdag(skel), retry = udag2pdagSpecial(skel)$pcObj, 
           relaxed = udag2pdagRelaxed(skel))
  }
  else {
    if (u2pd != "relaxed") 
      stop("Conservative PC can only be run with 'u2pd = relaxed'")
    tmp <- pc.cons.intern(skel, suffStat, indepTest, alpha, 
                          verbose = verbose, version.unf = c(2, 2))
    tripleList <- tmp$unfTripl
    if ( length(tripleList) >= 1 ) return(NULL)
    #print(tripleList)
    #browser()
    udag2pdagRelaxed(gInput = skel, verbose = verbose, unfVect = tripleList)
  }
}
