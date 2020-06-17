#Slightly modified version from the package deal.
#See http://cran.r-project.org/web/packages/deal/index.html for references.

jointprior.mod<-function (nw, N = NA, phiprior = "bottcher", timetrace = FALSE) 
{
    if (timetrace) {
        t1 <- proc.time()
        cat("[Jointprior ")
    }
    if (nw$nd > 0) {
        jointprob <- jointdisc(nw, timetrace = timetrace)
        minN <- min(2/jointprob)
        if (is.na(N)) 
            N <- minN
        if (N < minN) {
            cat("Warning: Your choice of imaginary sample size is very low\n")
            cat("We advice you to set the imaginary sample size to more than", 
                minN, "\n")
        }
        cat("Imaginary sample size:", N, "\n")
        jointalpha <- jointprob * N
        jointnu <- jointalpha
        jointrho <- jointalpha
    }
    else {
        jointnu <- N
        jointrho <- N
        jointalpha <- N
    }
    if (nw$nc > 0) {
        NN <- prod(dim(jointalpha))
        if (nw$nd > 0) {
            Dim <- dim(jointalpha)
            dparents <- nw$discrete
            lvek <- c()
            for (i in 1:NN) {
                cf <- findex(i, Dim, FALSE)
                label <- ""
                for (j in 1:ncol(cf)) {
                  label <- paste(label, nw$nodes[[dparents[j]]]$levelnames[cf[1, 
                    j]], sep = ":")
                }
                lvek <- c(lvek, label)
            }
        }
        jointmu <- matrix(NA, NN, nw$nc)
        jointsigma <- list()
        jointphi <- list()
        jcont <- jointcont.mod(nw, timetrace = timetrace)
        jointmu <- jcont$mu
        jointsigma <- jcont$sigma2
        dnames <- colnames(jointmu)
        for (i in 1:NN) {
            if (phiprior == "bottcher") {
                jointphi[[i]] <- jointsigma[[i]] * (jointnu[i] - 
                  1)
            }
            else {
                if (phiprior == "heckerman") {
                  jointphi[[i]] <- (jointrho[i] - 2)/(jointnu[i] + 
                    1) * jointnu[i] * jointsigma[[i]]
                }
                else stop("No such phiprior implemented")
            }
            colnames(jointmu) <- dnames
            colnames(jointsigma[[i]]) <- dnames
            rownames(jointsigma[[i]]) <- dnames
            colnames(jointphi[[i]]) <- dnames
            rownames(jointphi[[i]]) <- dnames
        }
        if (nw$nd > 0) {
            names(jointsigma) <- lvek
            names(jointphi) <- lvek
            rownames(jointmu) <- lvek
        }
    }
    else {
        jointphi <- NA
        jointmu <- NA
        jointsigma <- NA
    }
    if (timetrace) {
        t2 <- proc.time()
        cat((t2 - t1)[1], "]\n")
    }
    list(jointalpha = jointalpha, jointnu = jointnu, jointrho = jointrho, 
        jointmu = jointmu, jointsigma = jointsigma, jointphi = jointphi)
}
