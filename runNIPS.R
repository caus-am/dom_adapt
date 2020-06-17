### Main loop for algorithms

# For now only simulated data.
runNIPS <- function(algorithms = c("fs", "fs_jci_indep"), 
                    # List of algorithms we wish to evaluate. Default:
                    # - fs = standard feature selection
                    # - fs_jci_indep = algorithm described in the paper (fs + check is a provably separating set)
                    experiment_colname = "experiment", 
                    # Name of the column that contains the dataset indicator 
                    # (separating the different environments)
                    # For example, "geno" for original mouse challenge data
                    test = "logp", 
                    # Conditional independence test, default partial correlation.
                    # Some other possibilities in tests/ folder, not used in the NeurIPS paper.
                    p = 0.05, # Significance level for the partial correlation test.
                    threshold = 0, # If the score is above the threshold, it's an independence
                    dataDir=NULL,  # If !is.null(dataDir), read data from the files [dataDir]/[real|nan]-[i].csv
                    howmany=200, # How many randomly generated graphs we evaluate on
                    ### Options for simulator:
                    n=3, # number of observed variables
                    numInts=2, # number of intervention (domains)
                    NObs=1000, # number of data points in the observational case
                    NPerInt=1000, # number of data points in each interventional case
                    typeOfSimulator="int_NIPS2017", # which simulator to use for the data.
                    filter_require_path_IY=FALSE, 
                    IFactor=1.0,
                    ### General options again:
                    useMouseChallenge = FALSE,
                    useFsBruteforce=TRUE, #FALSE,
                    # Parameters for the Python scripts with random forests:
                    n_estimators_eval=1000, n_estimators_fs=1000, # size of random forests
                    currentDir="tmp", 
                    verbose = FALSE) {
  
  # Enabling Python on UvA server.
  # pythonCmdIntro <- paste0("source /opt/rh/python27/enable; cd ~/python; source bin/activate; cd -;")
  pythonCmdIntro <- ""
  evalCmd <- paste0(pythonCmdIntro, "python ../py_scripts/eval.py")
  
  if (is.null(dataDir)) {
    if (!useMouseChallenge) {
      dataDir <- paste0(currentDir,'/data/')
      system(paste("mkdir -p", dataDir))
    } else {
      dataDir <- '../data'
    }
  }
  
  losses <- foreach (i=1:howmany) %dopar% {
    if (verbose) cat("\n", format(Sys.time(), format="%H:%M:%S"), "***** Job", i, "started *****\n\n")

    # Initialize losses for FS and our method:
    fs_loss <- fs_jci_indep_loss <- 0
  
    # Where to store the results for the model i:
    mainDir <- paste0(currentDir, "/", i, "/")
    system(paste("mkdir -p", mainDir))

    # Full simulation data (no missing values)
    simFullFilename <- paste0(dataDir,'/real-',i,'.csv')
    # Simulation data that is actually used (with missing values of YBlind for domain IBlind)
    simFilename <- paste0(dataDir,'/nan-',i,'.csv')
    # Metadata file (number of variables, interventions, etc.)
    metaFilename <- paste0(dataDir,'/meta-',i,'.csv')
    
    if (!is.null(dataDir)) {
      # Assumes the following:
      # - (keep this) the regime variable comes after all system variables and is named [experiment_colname] (default "experiment"); intervention variables if any come after that;
      # - if no intervention variables are provided, use a diagonal design, with R=1 the observational dataset, R=2 corresponding to C1=1 etc;
      # - regimes are always numbered consecutively starting with 1.
      metadata <- as.list(read.csv(metaFilename, header=TRUE))
      sys_vars <- metadata$num_sys_vars
      num_vars <- sys_vars + metadata$num_ints_vars + 1
      sim <- list(M = list(fullG = array(0, c(num_vars, num_vars)),
                           G = array(0, c(sys_vars, sys_vars)),
                           YBlind=metadata$YBlind, IBlind=metadata$IBlind))
    } else if (!useMouseChallenge) {
      set.seed(i)
      sim <- simulateData.NIPS2017(n=n, numInts=numInts, NObs=NObs, NPerInt=NPerInt, typeOfSimulator=typeOfSimulator, IFactor=IFactor, filter_require_path_IY=filter_require_path_IY)
      if (!file.exists(simFilename)){
        write.csv(sim$D[[1]]$data_blinded, simFilename, row.names=FALSE)
      }
      if (!file.exists(simFullFilename)){
        write.csv(sim$D[[1]]$data, simFullFilename, row.names=FALSE)
      }
      num_vars <- dim(sim$M$fullG)[1]
      sys_vars <- dim(sim$M$G)[1]

      metadata <- data.frame(num_sys_vars=sys_vars, num_ints_vars=numInts,
                             YBlind=sim$M$YBlind, IBlind=sim$M$IBlind, IBlind_visible_value=0)
      write.csv(metadata, metaFilename, row.names=FALSE)
    } else {
      # Special case for mouse challenge (3 system variables, 2 intervention variables + observational case):
      sys_vars <- 3
      sim <- list(M = list(fullG=array(0, c(sys_vars+3, sys_vars+3)), 
                          G = array(0, c(sys_vars, sys_vars)), 
                          YBlind=3, IBlind=1))
    }
    regime_index <- sys_vars + 1
    
    # Binary encoding of sets used in ASP encoding.
    allSetsMinus_YBlindRIBlind <-  (2**(num_vars) - 1) - (2**(regime_index-1)) - (2**(sim$M$YBlind-1)) -  (2**(regime_index + sim$M$IBlind-1))
    
    # Run Feature selection and get top set.
    if ("fs" %in% algorithms || "fs_jci_indep" %in% algorithms) {
      fsFilename <- paste0(mainDir, "fs-",i,"-",sim$M$YBlind, "-", sim$M$IBlind, "-", allSetsMinus_YBlindRIBlind, ".csv")
      if (!file.exists(fsFilename)){
        if (useFsBruteforce) {
          fsRun <- paste0(pythonCmdIntro, "python ../py_scripts/fs_bruteforce.py") 
        } else {
          fsRun <- paste0(pythonCmdIntro, "python ../py_scripts/fs.py")
        }
        system(paste(fsRun, simFilename, sim$M$YBlind, sim$M$IBlind, allSetsMinus_YBlindRIBlind, mainDir, n_estimators_fs, sep = " "))
      }
      topFs <- as.list(read.csv(fsFilename,header=FALSE))[1]
    }
    
    ### Baseline method:
    if ("fs" %in% algorithms){  
      # Evaluate feature selection loss
      fsEvalFilename <- paste0(mainDir, "eval-",i,"-",sim$M$YBlind, "-", sim$M$IBlind, "-", topFs, ".csv")
      if (!file.exists(fsEvalFilename)){ 
        system(paste(evalCmd,simFullFilename, sim$M$YBlind, sim$M$IBlind, topFs, mainDir, n_estimators_eval, sep = " "))
      }
      fs_loss <- as.numeric(read.csv(fsEvalFilename,header=FALSE))
    }
    
    ### Our method:
    if ("fs_jci_indep" %in% algorithms) {
      jciIndepResult <- runJciIndep(i=i, simulatedData=sim, simFilename=simFilename, 
                                    experiment_colname = experiment_colname,
                                    threshold=threshold, test=test, p=p, mainDir=mainDir, 
                                    featureSelectionOutputFile=fsFilename, verbose=verbose)
      indeps<-jciIndepResult$indeps
      
      if (is.na(indeps)) {
        fs_jci_indep_loss <- NA
      } else {
        fsJciEvalFilename <- paste0(mainDir, "eval-",i,"-",sim$M$YBlind, "-", sim$M$IBlind, "-", indeps, ".csv")
        if (!file.exists(fsJciEvalFilename)){
          system(paste(evalCmd,simFullFilename, sim$M$YBlind, sim$M$IBlind, indeps, mainDir, n_estimators_eval, sep = " "))
        }
        fs_jci_indep_loss <- as.numeric(read.csv(fsJciEvalFilename,header=FALSE))
      }
    }
    
    if (verbose) {
      cat("\n", format(Sys.time(), format="%H:%M:%S"), "***** Job", i, "completed *****\n")
      cat("fs_loss =", fs_loss, "\n", "fs_jci_indep_loss =", fs_jci_indep_loss, ", score=", jciIndepResult$score,"\n")
    }
    
    list(fs_loss=fs_loss, fs_jci_indep_loss=fs_jci_indep_loss, indep_score=jciIndepResult$score)
  } # end foreach(i=1:howmany)
  
  if (verbose) cat("\n***** ALL JOBS COMPLETED *****\n")

  fs_losses <- sapply(losses, function(x) {x$fs_loss})
  fs_jci_indep_losses <- sapply(losses, function(x) {x$fs_jci_indep_loss})

  # Fallback to FS in case of NA
  fs_jci_nona_losses <- fs_jci_indep_losses
  fs_jci_nona_losses[which(is.na(fs_jci_indep_losses))] <- fs_losses[which(is.na(fs_jci_indep_losses))]
  
  # Feature selection losses for cases in which FS_JCI find something.
  fs_when_hej_nona_losses <- NA
  fs_when_hej_nona_losses[which(!is.na(fs_jci_indep_losses))] <- fs_losses[which(!is.na(fs_jci_indep_losses))]
  
  num_na_fs_jci <- sum(is.na(fs_jci_indep_losses))
  
  ldf <- data.frame(algorithm=c("Feature selection", "Our method"), loss=c(0,0), logloss=c(0,0), indep_score=c(0,0))
  for (l in losses) {
    ldf <- rbind (ldf, list(algorithm="Feature selection", loss=l$fs_loss, logloss=log(l$fs_loss), indep_score=0))
    ldf <- rbind (ldf, list(algorithm="Our method", loss=l$fs_jci_indep_loss, logloss=log(l$fs_jci_indep_loss), indep_score=l$indep_score))  }
  # Remove dummy rows.
  ldf <- ldf [c(-1,-2),]

  ldfFilename <- paste0(currentDir, "/ldf.csv")
  write.csv(ldf, ldfFilename, row.names=FALSE)

  #g<-ggplot(ldf, aes(algorithm, logloss)) + geom_boxplot(aes(algorithm, logloss), outlier.colour = "red", outlier.shape = 1) + xlab("Algorithm") + ylab("Mean squared loss (log-scale)")
  # ggplot(ldf, aes(algorithm, sqrt(loss))) + geom_boxplot(aes(algorithm, sqrt(loss)), outlier.colour = "red", outlier.shape = 1) + xlab("Method") + ylab("L2 loss") + coord_cartesian(ylim = c(0, 0.5)) 
  avgloss <- list(
             # mean of fs
             fs_loss=mean(fs_losses, na.rm=TRUE),
             # mean of fs_jci excluding NA
             nofallback_fs_jci_loss=mean(fs_jci_indep_losses, na.rm=TRUE),
             # mean of fs when fs_jci finds something != NA
             fs_when_hej_nona_losses = mean(fs_when_hej_nona_losses, na.rm=TRUE),
             # mean of fs_jci with fallback to FS when is NA
             fs_jci_loss=mean(fs_jci_nona_losses, na.rm=TRUE))
  list(avgloss = avgloss, all_losses=losses, ldf=ldf, num_na_fs_jci=num_na_fs_jci)
}
