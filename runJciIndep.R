## Run JCI and get all independences (default with HEJ).
runJciIndep <- function(i = 1, #seed offset for random generation
                    simulatedData=NULL, featureSelectionOutputFile=NULL,
                    encode="new_wmaxsat_acyclic_causes.pl",
                    threshold = 0, # If the score is above the threshold, it's an independence
                    experiment_colname = "experiment", # experiment_colname <- "geno" # for mouse challenge
                    mainDir=paste0("tmp/", i, "/"),
                    test="logp", p=0.05,
                    simFilename = paste0(mainDir,'/sim-', i, '-blinded.csv'),
                    verbose = FALSE) {
  system(paste("mkdir -p", mainDir))
  work_dir <- paste0(mainDir,"jci_indep/")
  system(paste("mkdir -p", work_dir))
  
  # Independence test data, if it exists, runJciIndep was already run on the same data.
  indepsRdata <- paste0(work_dir,"/indeps.Rdata")
  
  if (file.exists(indepsRdata)) {
    load(indepsRdata) # loads L
    #L <- list (indeps = cond_set)
    return (L)
  }
  
  if (is.null(simulatedData)) {
    # Simulate one model.
    set.seed(i)
    system(paste("mkdir -p", work_dir))
    sim <- simulateData.NIPS2017()
  } else {
    sim <- simulatedData
  }
  
  if (!file.exists(simFilename)){
    write.csv(sim$D[[1]]$data, simFilename, row.names=FALSE)
  }

  if (verbose) {
    fullG <- sign(sim$M$fullG + 1)
    sim$M$fullC <- dag2causes(fullG)
  }
  
  num_vars <- dim(sim$M$fullG)[1]
  sys_vars <- dim(sim$M$G)[1]
  
  IBlind <- sim$M$IBlind + sys_vars + 1
  
  cond_sets <- as.list(read.csv(file=featureSelectionOutputFile,header=FALSE))
  found <- FALSE
  scores <- c()
  scoreAccepted <- NA
  
  for (cond_set in cond_sets) {
      # Skip if cond_sets contains YBlind or IBlind
      bin_cond_set <- rev(dec.to.bin(cond_set, ndigits=num_vars))
      if (bin_cond_set[IBlind] == 1 || bin_cond_set[sim$M$YBlind] == 1 || bin_cond_set[sys_vars+1] == 1) {
          cat("runJciIndep skipping illegal feature set ", cond_set, "\n")
          next
      }
    
      if (verbose) {
        cat ("Test if intervention:", IBlind, " is independent of ", sim$M$YBlind, 
                        " given ", cond_set, "\n")
      }
      indepFilename <- paste0(work_dir, "/bckg_indep_", cond_set, ".txt")
      indepFile <- file(indepFilename)
      M <- rep(1,num_vars);
      M[c(which(bin_cond_set==1), sim$M$YBlind, IBlind)] <- 0
      M <- bin.to.dec(rev(M))
      writeLines( paste0(":- th(",IBlind, ",",sim$M$YBlind, ",",cond_set, ",0,",M,").\n", 
                         ":- th(",sim$M$YBlind, ",",IBlind,",",cond_set, ",0,",M,").\n",
                         ":- hh(",sim$M$YBlind, ",",IBlind, ",",cond_set, ",0,",M,").\n",
                         ":- tt(",sim$M$YBlind, ",",IBlind, ",",cond_set, ",0,",M,").\n",
                         ":- edge(",IBlind, ",", sim$M$YBlind,")."), indepFile)
  
      close(indepFile)

      result_indep <- run_experiment(filename=simFilename,
          n=num_vars, schedule=num_vars-2, solver="clingo", encode=encode, multipleMinima="baseline",
          experiment = experiment_colname, reference_dataset_label = 1, transform=NULL, test=test, p = p,
          background_knowledge_file=indepFilename, tmpDir=paste0(work_dir, "/tmp_indep_", cond_set), 
          verbose=verbose)
  
      depFilename <- paste0(work_dir, "/bckg_dep_", cond_set, ".txt")
      depFile <- file(depFilename)
      writeLines( paste0(":- not th(",IBlind, ",", sim$M$YBlind, ",",cond_set, ",0,",M,"),",
                         " not th(",sim$M$YBlind,  ",", IBlind, ",",cond_set, ",0,",M,"),",
                         " not hh(",sim$M$YBlind,  ",", IBlind, ",",cond_set, ",0,",M,"),",
                         " not tt(",sim$M$YBlind,  ",", IBlind, ",",cond_set, ",0,",M,").\n",  
                         ":- edge(",IBlind, ",", sim$M$YBlind,")."), depFile)
      close(depFile)

      result_dep <- run_experiment(filename=simFilename,
          n=num_vars, schedule=num_vars-2, solver="clingo", encode=encode, multipleMinima="baseline",
          experiment = experiment_colname, reference_dataset_label = 1, transform=NULL, test=test, p = p,
          background_knowledge_file=depFilename, tmpDir=paste0(work_dir, "/tmp_dep_", cond_set),
          verbose=verbose)
  
      if (result_dep$objective == -1) {
        #UNSAT 
        result_dep$objective = Inf
      }
      if (result_indep$objective == -1) {
        #UNSAT 
        result_indep$objective = Inf
      }
      score <- result_dep$objective - result_indep$objective
      
      scores <- c(scores, score)
      if (!is.na(score))
        if (score > threshold) {
          found <- TRUE
          scoreAccepted <- score
          break
        }
  }
  if (!found) cond_set = NA
  L <- list (indeps=cond_set, score=scoreAccepted)
  save(L, file=indepsRdata)
  L
}
