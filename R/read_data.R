run_experiment <- function(filename=NULL, # the filename of the csv with the data, NULL generates a simulated model.
                           indPath=NULL, # the filename of the in/dependence tests in the ASP format.
                           n = 3, # used only for the simulated model.
                           schedule=1, # Order of independence tests (max size of conditioning set)
                           N = 500, # used only for the simulated model.
                           confounder_proportion = 0.5, # used only for the simulated model.
                           numInts=3, # used only for the simulated model.
                           typeOfSimulator = "int", # used to decide which simulator to use.
                           solver="clingo", # You probably want this.
                           encode="new_wmaxsat_acyclic_causes.pl",
                           solver_conf="--quiet=1", # Configuration passed to clingo (quiet means less verbose possible).
                           multipleMinima="iterativeParallel", 
                           experiment = "experiment", # A column that represents the id of the experiment 
                           # In Sachs it's experiment, in MICE it's geno.
                           reference_dataset_label = 1, # The value the 'observational' dataset in the experiment column.
                           # In Sachs it's 1 (observational data, noICAM), in MICE it's 0 (wild type mice).
                           transform=NULL, # Transforms all data by log, if NULL doesn't change data
                           # In principle can be passed any function that takes a single argument.
                           test="logp", # Use logp test, since 'bayes' doesn't work well.
                           p = 0.05, # significance level in logp test.
                           only_vars=NULL, # if NULL takes all variables, otherwise specify which variables to use
                           # e.g. only_vars = c("raf", "mek") uses only the "raf" and "mek" columns.
                           only_ints=NULL, # if NULL takes all datasets, otherwise specify which datasets to use
                           # e.g. only_ints = c(1, 3) uses only the datasets that have values 1 and 3 in their 
                           # experiment column
                           evalDirectCauses=FALSE,
                           addBackgroundKnowledge=TRUE, # Add knowledge about regime and intervention variables.
                           background_knowledge_file=NULL, # Filename containing extra knowledge in ASP format, NULL if none.
                           postprocessingType="correct", 
		                       naiveVersion = FALSE,
		                       csv_sep = ",",
                           tmpDir="../tmp/log_sachs_data/", 
		                       verbose = FALSE) {
  system(paste("mkdir -p", tmpDir))

  if (is.null(indPath)) {
    M <- NULL
    if (is.null(filename)){
      simulationConfig <- list(n=n, topology="random", exconf='passive', N=N, restrict = "acyclic", confounder_proportion=confounder_proportion, numInts=numInts)
      MD <- simulateData(simulationConfig=simulationConfig, samples=NULL, model=NULL, indPath=NULL, returnData=TRUE, typeOfSimulator=typeOfSimulator) 
      filename <- paste(tmpDir, "simulated.csv", sep="")
      write.csv(x=MD$D[[1]]$data, file=filename, row.names=FALSE)
      plot_graph_as_dot.plotCauses(MD$M$fullC, filename = paste(tmpDir, "/true_causes_full.dot", sep=""), plotNotCauses=FALSE, verbose=F)
      tempM <- list(G=MD$M$fullG, Ge=MD$M$fullGe)
      plot_graph_as_dot.plot(tempM, filename = paste(tmpDir, "/true_full.dot", sep=""), plotNotCauses=FALSE, verbose=F)
      M <- MD$M 
      csv_sep <- ","
    }
    
    all_data <- preprocess_data(input=filename, background_knowledge_file=paste(tmpDir,"/bckg.txt", sep=""), 
                                currentDir=tmpDir,
                                p = p, alpha=1.15, test=test, schedule=schedule,
                                solver=solver,
                                only_vars_by_label=only_vars,
                                only_ints_by_label=only_ints,
                                reference_dataset_label = reference_dataset_label,
                                encode=encode,experiment = experiment,
                                csv_sep = csv_sep,
                                transform=transform,
                                M=M, addBackgroundKnowledge=addBackgroundKnowledge, 
                                postprocessingType=postprocessingType, verbose=verbose)
    samples <- all_data$data
    
    if (!is.null(background_knowledge_file)) {
      system(paste("cat", background_knowledge_file, ">>", all_data$bgk))
    }
    
    results <- pipeline(n=ncol(samples), schedule=schedule, encode=encode, samples=samples, 
                        background_knowledge_file=all_data$bgk,
                        indPath=all_data$indPath, tmpDir=tmpDir, model=M, evalDirectCauses=evalDirectCauses,
                        solver = solver, solver_conf=solver_conf, multipleMinima=multipleMinima, verbose=verbose)
  } else {
    results <- pipeline(schedule=schedule, encode=encode, indPath=indPath, 
                        background_knowledge_file=background_knowledge_file, tmpDir=tmpDir, evalDirectCauses=evalDirectCauses,
                        solver = solver, solver_conf=solver_conf, multipleMinima=multipleMinima, verbose=verbose)
    
  }
  #save(results, file=paste(tmpDir, "/results.Rdata", sep=""))
  return(results)
}

preprocess_data <- function(input="sachs2005/combined.csv",
                            currentDir = "../tmp/",
                            background_knowledge_file=paste(currentDir, "/bckg.txt", sep=""),
                            solver="clingo",
                            only_vars_by_label=NULL, # e.g. c('praf', 'pmek', 'p44.42')
                            only_ints_by_label=NULL, # labels in experiment.
                            reference_dataset_label = 1,
                            experiment = "experiment",
                            ints_labels = NULL,
                            csv_sep = "\t",
                            transform=NULL, 
                            p = 0.1,
                            alpha = 32,
                            test="bayes", schedule=1,
                            M = NULL,
			                      encode = "tom_acyclic_maxschedule_complete.pl",
                            addBackgroundKnowledge=TRUE, 
                            postprocessingType="correct",
                            simulateMissing=FALSE,
			                      verbose=FALSE) {
  all_data <- read.csv(input, header = TRUE, sep=csv_sep, quote = "\"",
                       dec = ".", fill = TRUE)
  
#   if (simulateMissing) {
#     all_data[which(all_data[[experiment]] == 2), 1] <- NA
#     all_data[which(all_data[[experiment]] == 3), 2] <- NA
#   }
   
  bckgFile <- file(background_knowledge_file, "w")
  
  if (!is.null(only_ints_by_label)) {
    all_data <- all_data[all_data[[experiment]] %in% only_ints_by_label, ]
  } else {
    only_ints_by_label <- sort(unique(all_data[[experiment]]))
  }
  
  if (!is.null(only_vars_by_label)) {
    # Note: this probably gives unexpected results if the data file included context variables
    all_data <- all_data[, c(only_vars_by_label, experiment)] 
  }
  
  num_vars_total <- dim(all_data)[2]
  standard_vars_col <- match(experiment,names(all_data)) - 1 # Assume all system variables are to the left of experiment variable
  
  if (!is.null(transform)){
    if (identical(transform, log) && any(all_data [,1:standard_vars_col] < 0)) {
      stop("Log transform on at least one negative value.")
    }
    all_data <- cbind(transform(all_data [,1:standard_vars_col]), all_data[,(standard_vars_col+1)])
  }

  if (num_vars_total == standard_vars_col + 1) {
    # Construct diagonal design
    count_intervention_vars <- 0
    if (length(only_ints_by_label) > 1) {
      for (i in only_ints_by_label) {
        if (i == reference_dataset_label) next
        if (!is.null(ints_labels)){
          int_name <- paste(ints_labels[i], sep="")
        } else{
          int_name <- paste("i", i, sep="")
        }
        all_data[int_name] <- 0
        all_data[all_data[[experiment]]==i, int_name] <- 1
        count_intervention_vars <- count_intervention_vars + 1
      }
    }
  } else {
    count_intervention_vars <- num_vars_total - standard_vars_col - 1
    cat("Total vars:", num_vars_total, "\n")
    cat("Sys   vars:", standard_vars_col, "\n")
    cat("Int   vars:", count_intervention_vars, "\n")
  }
  
  exp_index <- standard_vars_col + 1
  int_indices <- (standard_vars_col +  2 ) : (standard_vars_col + 1 + count_intervention_vars)
  standard_vars_indices <- 1:standard_vars_col
  n <- standard_vars_col + 1 + count_intervention_vars
  
  if (addBackgroundKnowledge && length(only_ints_by_label) > 1) {
    # Experiment causes all ints and nothing causes experiment:
    # Nothing causes intervention variables (except experiment).
    if (solver != "clingo") { 
      close(bckgFile)
      stop("Cannot work with other solvers.")
    } else {
      ### CLINGO Version.
      
      # Test: Standard vars || Standard vars | any* except (all indicators + R)
      # Test: Indicator vars || standard vars | any but R
      # Test: Regime var || Standard var | any but all indicators
      
      cat("vars(1..", standard_vars_col,"). ints(", min(int_indices), "..", max(int_indices),"). ", 
          "regime(", exp_index,").\n", sep="", file=bckgFile)
      
      if (grepl('new_', encode)){
         # Direct causal relation method 
         cat(":- not edge(R, I), regime(R), ints(I).\n", sep="", file=bckgFile) 
         cat(":- conf(R, A), regime(R), node(A).\n", sep="", file=bckgFile)
	       cat(":- conf(A, R), regime(R), node(A).\n", sep="", file=bckgFile)
         cat(":- edge(R, X), regime(R), vars(X).\n", sep="", file=bckgFile)
	       cat(":- edge(A, R), regime(R), node(A).\n", sep="", file=bckgFile)
      	 cat(":- conf(I, A), ints(I), node(A).\n", sep="", file=bckgFile)
      	 cat(":- conf(A, I), ints(I), node(A).\n", sep="", file=bckgFile)
         cat(":- edge(X, I), vars(X), ints(I).\n", sep="", file=bckgFile)
         cat(":- edge(I, J), ints(I), ints(J).\n", sep="", file=bckgFile)
         close(bckgFile)
      } else {     

      # Nothing causes R, I_i is caused only by R.
      cat(":- causes(X, I), vars(X), ints(I).\n",
          ":- causes(I, J), ints(I), ints(J).\n",
          ":- causes(X, R), regime(R), node(X).\n",
          ":- not causes(R,Y), regime(R), ints(Y).\n",
          sep="", file=bckgFile)

      cat("\n% Extra rules: \n", sep="", file=bckgFile)
      
      # 2. If R causes a var, then it is through one of the I_i.
      cat("{ causes(I,X): ints(I) } :- causes(R, X), vars(X), regime(R). \n", sep="", file=bckgFile)
      # 5. dep(R,X) => R --> X (no confounders).
      cat(":- dep(R, X, 0), not causes(R, X), regime(R), vars(X) .\n", sep="", file=bckgFile)
      # Possibly redundant:
      # indep(X,Y) => R -/-> X or R -/-> Y
      cat(":- indep(X, Y, 0), causes(R, X), causes(R, Y), regime(R), R!=X, R!=Y, X!=Y. \n", sep="", file=bckgFile)
      
      # (extra) Indicator var || Indicator var | R = indepedendent 
      #cat(":- dep(I, J, U), U==2**(R-1), regime(R), ints(I), ints(J), I != J.\n", sep="", file=bckgFile) 
      
      #cat("\n% Regime var x Indicator var \n", sep="", file=bckgFile)
      # (extra) Regime var || Indicator var | any = dependent (R is a direct cause of I_i)
      #cat(":- indep(R,J,Z), set(Z), regime(R), not ismember(Z,R), not ismember(Z,J), ints(J).\n", sep="", file=bckgFile)       
      
      ## 6. dep(I_i, X, {I_j} j=/=i) => I_i --> X  (no confounders) 
      #cat(":- dep(X, I, SUMOTHERS), not causes(I, X), vars(X), ints(I),",
      #    "SUMOTHERS == 2**(n) - 2**(I-1) - 2**(", n+1 ,"). \n", sep="", file=bckgFile) 
      cat(":- dep(X, I, R), not causes(I, X), vars(X), ints(I), regime(R).\n",sep="", file=bckgFile)
      
      cat("\n% 1. Regime var x Indicator var \n", sep="", file=bckgFile)
      # Regime var || Indicator var | any = dependent (R is a direct cause of I_i)
      cat(":- indep(R,J,Z), set(Z), regime(R), not ismember(Z,R), not ismember(Z,J), ints(J).\n", sep="", file=bckgFile)      

      # 3. Regime var || Standard var | any plus all indicators = independent
      cat("\n% Regime var x Standard var | any plus all indicators\n", sep="", file=bckgFile)
      cat(":- dep(X, R, U), set(U), U & SUMINTS == SUMINTS, set(SUMINTS), SUMINTS == 2**(n) - 2**(R), ",
          " not ismember(U, R), not ismember(U, X), vars(X), regime(R). \n", sep="", file=bckgFile)  

      cat("\n% Cannot have minDep with R. \n", sep="", file=bckgFile)
      # 7, 8. Cannot have minDep with I_i or R.
      # :- indep(X, Y, Z), dep(X,Y, Z u R).
      cat(":- indep(X,Y,W), dep(X,Y,U), U==W+2**(R-1), regime(R), node(X), node(Y), set(W), set(U),",
          " not ismember(W, R), not ismember(W, X), not ismember(W, Y),", 
          " not ismember(U, X), not ismember(U, Y),",
          " X != Y, X != R, Y != R.\n", sep="", file=bckgFile)
      
      cat("\n% Cannot have minDep with I_i.\n", sep="", file=bckgFile)
      # :- indep(X, Y, Z), dep(X,Y, Z u I_i).
      cat(":- indep(X,Y,W), dep(X,Y,U), U==W+2**(I-1), set(W), set(U),",
          " not ismember(W, I), not ismember(W, X), not ismember(W, Y),", 
          " not ismember(U, X), not ismember(U, Y), ", 
          " ints(I), node(X), node(Y), X != Y, X != I, Y != I.\n", sep="", file=bckgFile)
      close(bckgFile)
      }
    }
  }

  
  test_data <- list()
  test_data$X <- all_data
  
  if (test == "bayes") {
    test_data$p_threshold<-p # prior probability of ind.
    test_function<-test.bayes 
    test_data$alpha<-alpha # eq. sample size for the prior
    test_data$N<-nrow(all_data)
    test_data$discrete<-rep(0,n)
    test_data$discrete[exp_index] <- 1
    test_data$discrete[int_indices] <- 1
    test_data$discrete <- as.logical(test_data$discrete)
    df<-data.frame(all_data) 
    df[, test_data$discrete] <- lapply(df[, test_data$discrete], as.factor)
    test_data$df <- df
  } else if (test == "logp") {
    #test_data$Cx<-cov(all_data, use="pairwise.complete.obs")
    test_data$Cx<-cov(all_data)
    test_data$N<-nrow(all_data)
    test_data$p_threshold<-p # significance level.
    if (any(is.na(all_data))){
      test_function<-test.logp.missing
    } else {
      test_function<-test.logp
    }
  } else if (test == "oracle") {
    M <- list(G=M$fullG, Ge=M$fullGe)
    test_data$M<-M #should take out the Js here  
    test_data$N<-Inf
    test_data$p_threshold<-0.5 # independent vars have p-value 1, dependent 0.
    test_function<-test.oracle
  }
  
  data <- list()
  data$n <- n
  data$exp_index <- exp_index
  data$standard_vars_indices <- standard_vars_indices
  data$int_indices <- int_indices
  data$data <- all_data
  data$bgk <- background_knowledge_file
  
  tic()
# tested_independences <- conductIndepTests_mod.loop(test_function, test_data, maxcset=schedule, n=n, 
#                                                      standard_vars_indices, int_indices, exp_index)
  tested_independences_additional_bckg <- conductIndepTests_mod_naive.loop(test_function, test_data, maxcset=schedule, n=n, 
                                                       standard_vars_indices, int_indices, exp_index, 
                                                       postprocessingType=postprocessingType, verbose=verbose)
  tested_independences <- tested_independences_additional_bckg$tested_independences
  additional_bckg <- tested_independences_additional_bckg$additional_bckg
  test_time <- toc()
  
  indPath <- paste(currentDir, "/results.ind", sep='')
  writeIndepsToFile(n=n, tested_independences=tested_independences, write_function=write_constraint, 
                    write_data=list(weight="log"), indFilenameFullPath = indPath)   
  
  indFile <- file(indPath, "a")
  cat(additional_bckg, file=indFile)
  close(indFile)  
  writeAspSets.old(n=n, aspSetsFullPath = background_knowledge_file)
  build_full(n=n, aspSetsFullPath = background_knowledge_file,tested_independences=tested_independences)

  data$tested_independences <- tested_independences
  data$indPath <- indPath
  return(data)
}

conductIndepTests_mod.loop <- function(test_function, test_data, maxcset=Inf, n, 
                                       standard_vars_indices, int_indices, exp_index) {
  tested_independences <- list()
  
  # TODO: check if they are deterministic at the test side.
  
  cat("Test: Standard vars || Standard vars | any*\n")
  # Standard vars || Standard vars | any*
  for (csetsize in index(0,maxcset) ) {       
    for ( i in standard_vars_indices) {
      jindices <- standard_vars_indices[which(standard_vars_indices > i)]
      tested_independences_j <- foreach (j = jindices) %do% {
        conductIndepTests.parallel.loop(test_function, test_data, n, i, j, csetsize, conditioning_vars=NULL, except_all_vars_in_cond_set=c(int_indices, exp_index))
      } #for j
      for (test_result_list in tested_independences_j) {
        for (test_result in test_result_list) {
          test_result$jset <- 0
          tested_independences[[length(tested_independences) + 1]] <- test_result
        }        
      }
    } # for i
  } # for csetsize
  
  cat("Test: Indicator vars || Standard vars | any but R\n")
  # Indicator vars || standard vars | any but R
  for (csetsize in index(0,maxcset) ) {       
    for ( i in int_indices) {
      tested_independences_j <- foreach (j = standard_vars_indices) %do% {
        conductIndepTests.parallel.loop(test_function, test_data, n, j, i, csetsize, conditioning_vars=c(standard_vars_indices,int_indices))
      } #for j
      for (test_result_list in tested_independences_j) {
        for (test_result in test_result_list) {
          test_result$jset <- 0
          tested_independences[[length(tested_independences) + 1]] <- test_result
        }        
      }
    } # for i
  } # for csetsize
  
  cat("Test: Indicator vars || Indicator vars | any but R\n")
  # Indicator vars || indicator vars | any but R
  for (csetsize in index(0,maxcset) ) {       
    for ( i in int_indices) {
      jindices <- int_indices[which(int_indices > i)]
      tested_independences_j <- foreach (j = jindices) %do% {
        conductIndepTests.parallel.loop(test_function, test_data, n, j, i, csetsize, conditioning_vars=c(standard_vars_indices,int_indices))
      } #for j
      for (test_result_list in tested_independences_j) {
        for (test_result in test_result_list) {
          test_result$jset <- 0
          tested_independences[[length(tested_independences) + 1]] <- test_result
        }        
      }
    } # for i
  } # for csetsize
  
  cat("Test: Regime var || Standard var | any but all indicators\n")
  # Regime var || Standard var | any but all indicators is ok (the all indicators case is taken care before)
  for (csetsize in index(0,maxcset) ) {       
    tested_independences_j <- foreach (j = standard_vars_indices) %do% {
      conductIndepTests.parallel.loop(test_function, test_data, n, j, exp_index, csetsize, conditioning_vars=c(standard_vars_indices,int_indices), except_all_vars_in_cond_set=int_indices)
    } #for j
    for (test_result_list in tested_independences_j) {
      for (test_result in test_result_list) {
        test_result$jset <- 0
        tested_independences[[length(tested_independences) + 1]] <- test_result
      }        
    }
  } # for csetsize
  
  tested_independences
}  

conductIndepTests_mod_naive.loop <- function(test_function, test_data, maxcset=Inf, n, 
                                       standard_vars_indices, int_indices, exp_index, 
                                       postprocessingType="correct", verbose=FALSE) {
  tested_independences <- list()
  
  experimental_design_matrix <- unique(test_data$X[, c(exp_index, int_indices)])
  DL <- find_deterministic_relations(experimental_design_matrix, enforceMinimality = FALSE, exp_index=exp_index)
  
  all_vars <- c(standard_vars_indices,int_indices,exp_index)
  
  additional_bckg <- ""
  
  # Any x Any | Any
  for (csetsize in index(0,maxcset) ) {       
    for ( i in all_vars) {
      jindices <- all_vars[which(all_vars > i)]
      tested_independences_j <- foreach (j = jindices) %do% {
        conductIndepTests.parallel.loop(test_function, test_data, n, i, j, csetsize, conditioning_vars=NULL)
      } #for j
      for (test_result_list in tested_independences_j) {
        for (test_result in test_result_list) {
          test_result$jset <- 0
          if (length(test_result$C) > 0)  {
            mapped_variables <- test_result$C - exp_index + 1 
            mapped_variables <- mapped_variables[mapped_variables>0]
            determined <- unlist(getDet(mapped_variables, DL))
            
             if (length(determined)!= length(mapped_variables) || !all(determined == mapped_variables)) {
              if (verbose) cat("\nSkipping test result: ", test_result$vars, "|", test_result$C, "\n")
              
              remapped_determined <- determined + exp_index - 1
              remapped_determined <- c(remapped_determined, test_result$C)
              remapped_determined <- unique(remapped_determined)
              
              if ((postprocessingType == "wrong") && all(test_result$vars %in% remapped_determined)){
                test_result$independent <- TRUE
              } else {
                remapped_determined_cset <- rep(0, n)
                remapped_determined_cset[remapped_determined] <- 1
                remapped_determined_cset<-bin.to.dec(rev(remapped_determined_cset))
                
                if (postprocessingType == "wrong") {
                  additional_bckg <- paste(additional_bckg, 
                                         ":- indep(", test_result$vars[1], ",",  test_result$vars[2], ",",  test_result$cset, 
                                         "), dep(", test_result$vars[1], ",",  test_result$vars[2], ",",  remapped_determined_cset, ").\n", sep="")
                  additional_bckg <- paste(additional_bckg, 
                                           ":- dep(", test_result$vars[1], ",",  test_result$vars[2], ",",  test_result$cset, 
                                           "), indep(", test_result$vars[1], ",",  test_result$vars[2], ",",  remapped_determined_cset, ").\n", sep="")
                  
                }
                if (postprocessingType != "oracle") {
                  test_result$independent <- NA
                }
              }
            }
          }
          tested_independences[[length(tested_independences) + 1]] <- test_result
        }        
      }
    } # for i
  } # for csetsize
  list(tested_independences=tested_independences, additional_bckg=additional_bckg)
}  
