#Slightly modified version from the package deal.
#See http://cran.r-project.org/web/packages/deal/index.html for references.

networkfamily.mod <-function (data, nw = network(data), prior = jointprior(nw), trylist = vector("list",
    size(nw)), timetrace = TRUE)                                                     
{                                                                                    
    if (timetrace) {                                                                 
        t1 <- proc.time()                                                            
        cat("[networkfamily ")                                                       
    }                                                                                
    nw <- learn(nw, data, prior, trylist = trylist)                                  
    trylist <- nw$trylist                                                            
    nw <- nw$nw                                                                      
    ndiscrete <- nw$nd                                                               
    ncontinuous <- nw$nc                                                             
    cat("Creating all (", numbermixed(ndiscrete, ncontinuous),                       
        " minus restrictions) networks with ", ndiscrete, " discrete and ",          
        ncontinuous, " continuous nodes\n", sep = "")                                
    nwl <- list()                                                                    
    n <- ndiscrete + ncontinuous                                                     
    nwl <- list(nw)                                                                  
    for (node in 2:n) {                                                              
        for (idx in 1:length(nwl)) {                                                 
            nws <- addarrows(nwl[[idx]], node, data, prior, trylist = trylist)       
            trylist <- nws$trylist                                                   
            nwl <- c(nwl, nws$nw)                                                    
        }                                                                            
    }
    cat("Created", length(nwl), "networks, ")
    if (ndiscrete > 2 | ncontinuous > 2) {
        cat("removing cycles...\n")
        nwlres <- nwl[!unlist(lapply(nwl, cycletest))]
        cat(length(nwl) - length(nwlres), "cycles removed, ending up with",
            length(nwlres), "networks\n")
    }
    else nwlres <- nwl
    class(nwlres) <- "networkfamily"
    if (timetrace) {
        t2 <- proc.time()
        cat((t2 - t1)[1], "]\n")
    }
    list(nw = nwlres, trylist = trylist)
}
