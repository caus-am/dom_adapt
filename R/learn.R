learn <- function(testConfig, solverConfig, currentDir, graph_filename_template, indPath, background_knowledge_file, 
                  tested_independences, MD, evalDirectCauses=FALSE, verbose=0){ 
  # Function that learns the model.
  # Returns the learned model L.
  # available_solvers <- c(learn.asp.available_solvers, learn.pcalg.fci.available_solvers)
  if (solverConfig$solver %in% learn.pcalg.fci.available_solvers) {
    # TODO: Add fixedGaps = fixedGaps, fixedEdges = fixedEdges
    L <- learn.pcalg.fci(n=testConfig$n, solver=solverConfig$solver, alpha=testConfig$p, schedule=testConfig$schedule, 
                         tested_independences=tested_independences);
  } else if (solverConfig$solver %in% learn.asp.available_solvers) {   
    # Write independence tests results to file.
    if (is.null(indPath)) {
      indPath <- paste(currentDir, "/", graph_filename_template, '.ind', sep='')
      writeIndepsToFile(n=testConfig$n, tested_independences=tested_independences, write_function=write_constraint, 
                        write_data=list(weight=testConfig$weight), indFilenameFullPath = indPath)    
      if (length(solverConfig$intervened_variables) > 0) {
        # Check causal oracle and write causes and not causes relations.
        write_causal_constraint(MD$M$C, append=TRUE, ints = solverConfig$intervened_variables, indFilenameFullPath=indPath)
      }
    }
    L <- learn.asp(testConfig=testConfig, solverConfig=solverConfig, graph_filename_template=graph_filename_template,
                   currentDir=currentDir, indPath=indPath, background_knowledge_file=background_knowledge_file, 
                   tested_independences=tested_independences, evalDirectCauses=evalDirectCauses, verbose=verbose)
  } else if (solverConfig$solver == "bnlearn"){
    L <- learn.bnlearn(D=MD$D, n=testConfig$n, p_threshold=1) 
  } else{
    stop("Solver ", solver, " not found.")
  }
  
  L
}
